#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import joblib
from scipy.special import expit as sigmoid


def compute_pseudocount(series_list, scale_factor=0.1, min_pc=1e-6):
    mins = []
    for s in series_list:
        arr = np.array(s)
        pos = arr[arr > 0]
        if len(pos) > 0:
            mins.append(pos.min())
    if len(mins) == 0:
        base = min_pc
    else:
        base = max(min(mins) * scale_factor, min_pc)
    return float(base)


def build_feature_df(male, female, pseudocount):
    df = pd.DataFrame({
        "male_expr": male.astype(float),
        "female_expr": female.astype(float),
    })

    df["log2FC"] = np.log2((df["male_expr"] + pseudocount) / (df["female_expr"] + pseudocount))
    df["abs_log2FC"] = df["log2FC"].abs()
    df["ratio"] = (df["male_expr"] + pseudocount) / (df["female_expr"] + pseudocount)

    df["male_only"] = ((df["male_expr"] > 0.05) & (df["female_expr"] < 0.01)).astype(int)
    df["female_only"] = ((df["female_expr"] > 0.05) & (df["male_expr"] < 0.01)).astype(int)

    df["in_SDR"] = 1

    return df


def match_features(X, required_features):
    X2 = pd.DataFrame(index=X.index)
    for f in required_features:
        if f in X.columns:
            X2[f] = X[f]
        else:
            X2[f] = 0.0
    return X2


def compute_expression_confidence(expr_mean, median_expr, slope=0.3, center_scale=0.5):
    center = median_expr * center_scale
    return sigmoid(slope * (expr_mean - center))


def main():
    parser = argparse.ArgumentParser(description="Sex determination gene prediction for T.pilosa")
    parser.add_argument("--input", default="tpi.txt", help="Input expression file (tab-separated)")
    parser.add_argument("--model", default="model_final.joblib", help="Trained model files")
    parser.add_argument("--output", default="tpi_predictions_final.csv", help="Outputs predictions to this file")
    parser.add_argument("--expr_slope", type=float, default=0.3, help="Expression level confidence level sigmoid slope (smaller value = gentler filtering)")
    parser.add_argument("--expr_center_scale", type=float, default=0.5, help="Expression level center point scaling factor (smaller value = retain more low-expression genes)")
    args = parser.parse_args()

    obj = joblib.load(args.model)
    pipeline = obj["pipeline"]
    feature_list = obj["features"]

    df = pd.read_csv(args.input, sep="\t", dtype={"gene_id": str})
    needed = ["gene_id", "M1", "M2", "F1", "F2"]
    for c in needed:
        if c not in df.columns:
            raise ValueError(f"输入文件缺少必需列: {c}")

    pseudocount = compute_pseudocount([df['M1'], df['M2'], df['F1'], df['F2']])

    expr_mean_all = df[["M1", "M2", "F1", "F2"]].mean(axis=1)
    median_expr = np.median(expr_mean_all.values)

    feats1 = build_feature_df(df["M1"].values, df["F1"].values, pseudocount)
    feats1 = match_features(feats1, feature_list)
    prob1 = pipeline.predict_proba(feats1)[:, 1]

    feats2 = build_feature_df(df["M2"].values, df["F2"].values, pseudocount)
    feats2 = match_features(feats2, feature_list)
    prob2 = pipeline.predict_proba(feats2)[:, 1]

    combined_score = np.maximum(prob1, prob2)

    expr_conf = compute_expression_confidence(
        expr_mean_all, median_expr, 
        slope=args.expr_slope, 
        center_scale=args.expr_center_scale
    )

    final_score = combined_score * expr_conf

    out = pd.DataFrame({
        "gene_id": df["gene_id"],
        "M1": df["M1"],
        "M2": df["M2"], 
        "F1": df["F1"],
        "F2": df["F2"],
        "model_score_M1F1": prob1,
        "model_score_M2F2": prob2,
        "combined_score": combined_score,
        "expr_mean": expr_mean_all,
        "expr_confidence": expr_conf,
        "final_score": final_score,
    })

    for col in feats1.columns:
        out[f"M1F1_{col}"] = feats1[col]
    for col in feats2.columns:
        out[f"M2F2_{col}"] = feats2[col]

    out = out.sort_values("final_score", ascending=False)
    out["rank"] = range(1, len(out) + 1)

    out.to_csv(args.output, index=False)
    print(f"Prediction complete! Results saved to: {args.output}")

if __name__ == "__main__":
    main()
