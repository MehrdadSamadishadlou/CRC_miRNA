# Model hyperparameters present in the supplied code

## First layer: CRC versus non-cancer

```python
xgb.XGBClassifier(
    objective="binary:logistic",
    random_state=1
)
```

All other XGBoost parameters use the defaults of the installed XGBoost version.

Additional settings:

- Five-fold cross-validation.
- CRC is the positive class.
- Classification threshold: 0.50.
- Bootstrap confidence intervals: 2,000 resamples.
- Bootstrap random seed: 42.
- SHAP background sample: up to 200 training samples.
- SHAP sampling seed: 42.

## Second layer: low-risk versus high-risk CRC

Feature ranking and subset selection:

```python
xgb.XGBClassifier(
    objective="binary:logistic",
    random_state=1
)
```

Final model:

```python
xgb.XGBClassifier(
    n_estimators=1000,
    max_depth=7
)
```

Class-dependent sample weights:

```python
{
    0: 0.9114252061248527,
    1: 9.082451253481894
}
```

Additional settings:

- Low risk, stages 0–2: class 0.
- High risk, stages 3–4: class 1.
- Candidate subset selection: five-fold cross-validation with macro F1.
- Final classification threshold: 0.50.
- Bootstrap confidence intervals: 2,000 resamples.
- Bootstrap random seed: 42.
- SHAP background sample: up to 200 training samples.
- SHAP sampling seed: 42.

## Reproducibility warning

The final risk model in the supplied code does not set `random_state`. This was preserved rather than silently changing the model. The exact XGBoost version and execution environment therefore need to be recorded, and the author should decide whether adding a fixed seed is scientifically acceptable before release.
