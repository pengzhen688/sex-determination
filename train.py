#!/usr/bin/env python3
"""
Train script for plant sex-determining gene prediction
- Uses BalancedRandomForestClassifier from imblearn to handle extreme imbalance
- Compatible with predict.py
"""

import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score
from joblib import dump

def load_data(path):
    df = pd.read_csv(path, sep="\t")
    df['male_expr'] = pd.to_numeric(df['male_expr'], errors='coerce').fillna(0)
    df['female_expr'] = pd.to_numeric(df['female_expr'], errors='coerce').fillna(0)
    df['in_SDR'] = df['in_SDR'].astype(int)
    df['label'] = df['label'].astype(int)
    return df

def featurize(df):
    eps = 1e-6
    df = df.copy()
    df["log2FC"] = np.log2((df["male_expr"]+eps)/(df["female_expr"]+eps))
    df["abs_log2FC"] = df["log2FC"].abs()
    df["ratio"] = (df["male_expr"]+eps)/(df["female_expr"]+eps)
    df["male_only"] = ((df["male_expr"] > 0.05) & (df["female_expr"] < 0.01)).astype(int)
    df["female_only"] = ((df["female_expr"] > 0.05) & (df["male_expr"] < 0.01)).astype(int)
    feature_cols = ["male_expr","female_expr","log2FC","abs_log2FC","ratio","male_only","female_only","in_SDR"]
    X = df[feature_cols].astype(float)
    y = df["label"].values
    meta = df[["species","gene_id","annotation","label"]]
    return X, y, meta, feature_cols

def precision_at_k(y_true, scores, k):
    idx = np.argsort(scores)[-k:]
    return precision_score(y_true[idx], np.ones(k))

def evaluate_model(X, y, meta, model, prefix=""):
    scores = model.predict_proba(X)[:,1]
    roc = roc_auc_score(y, scores)
    ap = average_precision_score(y, scores)
    k = max(1, sum(y))
    p_at_k = precision_at_k(y, scores, k)
    print(f"{prefix}ROC AUC = {roc:.4f}, PR AUC = {ap:.4f}, Precision@k(k={k}) = {p_at_k:.4f}")
    df = meta.copy()
    df["score"] = scores
    df["rank"] = df["score"].rank(ascending=False, method="min").astype(int)
    tp = df[df["label"]==1].sort_values("rank")
    print(f"{prefix}Detected {tp.shape[0]} true positives:")
    print(tp[["gene_id","species","score","rank","annotation"]].to_string(index=False))

def main(args):
    df = load_data(args.input)
    X, y, meta, feature_cols = featurize(df)

    # BalancedRandomForest pipeline
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", BalancedRandomForestClassifier(
            n_estimators=500,
            max_depth=3,
            random_state=42
        ))
    ])

    # Stratified 3-fold CV
    print("\n=== Stratified 3-fold CV ===")
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
    cv_scores = []
    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y),1):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        model = pipe.fit(X_train, y_train)
        scores = model.predict_proba(X_test)[:,1]
        roc = roc_auc_score(y_test, scores)
        ap = average_precision_score(y_test, scores)
        cv_scores.append((roc, ap))
        print(f"Fold {fold}: ROC={roc:.4f}, AP={ap:.4f}")
    print("Mean ROC =", np.mean([s[0] for s in cv_scores]))
    print("Mean AP  =", np.mean([s[1] for s in cv_scores]))

    # Train on full data
    model_final = pipe.fit(X, y)

    # Train-set evaluation
    print("\n=== Train-set Evaluation (full data) ===")
    evaluate_model(X, y, meta, model_final, prefix="Train ")

    # LOSO evaluation
    print("\n=== Leave-One-Species-Out ===")
    for sp in df["species"].unique():
        mask = (df["species"] != sp)
        X_train, X_test = X[mask], X[~mask]
        y_train, y_test = y[mask], y[~mask]
        model = pipe.fit(X_train, y_train)
        if len(np.unique(y_test)) < 2:
            print(f"[{sp}] only 1 class in test, skipping metrics.")
            continue
        scores = model.predict_proba(X_test)[:,1]
        roc = roc_auc_score(y_test, scores)
        ap  = average_precision_score(y_test, scores)
        print(f"[{sp}] ROC={roc:.4f}, PR AUC={ap:.4f}")

    # Save model for predict.py
    dump({"pipeline": model_final, "features": feature_cols}, "model_final.joblib")
    print("\nSaved model_final.joblib (pipeline + features)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="combined_training.tsv")
    args = parser.parse_args()
    main(args)

