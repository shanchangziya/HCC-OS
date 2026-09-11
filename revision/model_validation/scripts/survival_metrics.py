"""Plot-free survival utilities, matching timeROC marginal IPCW point estimates.

Cases have T < horizon, controls have T > horizon. Individuals at exactly the
horizon are excluded as in timeROC. Reverse Kaplan-Meier weights use G(T-).
"""
import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored
from sksurv.nonparametric import kaplan_meier_estimator


def censor_weights(time, event):
    t, g = kaplan_meier_estimator(event.astype(bool), time, reverse=True)
    before = np.r_[1., g][np.searchsorted(t, time, side='left')]
    return np.divide(1., before, out=np.full_like(before, np.nan), where=before>0)


def roc_auc(time, event, score, horizon, return_curve=False):
    time, event, score = np.asarray(time), np.asarray(event, bool), np.asarray(score)
    cases = (time < horizon) & event
    controls = time > horizon
    if cases.sum() == 0 or controls.sum() == 0:
        return (np.nan, None) if return_curve else np.nan
    weights = censor_weights(time, event)
    wc = np.where(cases, weights, 0.)
    wn = controls.astype(float)
    order = np.argsort(-score, kind='stable')
    # Curves after each complete block of tied markers.
    last = np.r_[score[order][:-1] != score[order][1:], True]
    tp = np.r_[0., np.cumsum(wc[order])[last] / wc.sum()]
    fp = np.r_[0., np.cumsum(wn[order])[last] / wn.sum()]
    auc = np.trapezoid(tp, fp)
    if not return_curve:
        return auc
    threshold = np.r_[np.inf, score[order][last]]
    return auc, pd.DataFrame({'threshold': threshold, 'fpr': fp, 'tpr': tp})


def cindex(time, event, score):
    return concordance_index_censored(np.asarray(event, bool), np.asarray(time), np.asarray(score))[0]


def bootstrap_metrics(time, event, scores, times=(365., 1095., 1825.), B=1000, seed=20260905):
    """Paired sample-row bootstrap; fixed model scores, no retraining in these CIs."""
    time, event = np.asarray(time), np.asarray(event, bool)
    names = list(scores)
    matrix = np.column_stack([scores[n] for n in names])
    n = len(time)
    rng = np.random.default_rng(seed)
    # Precompute comparable pairs and all marker comparisons. Sample weights
    # (multinomial counts) exactly represent resampling complete sample rows.
    comparable = event[:, None] & ((time[:, None] < time[None, :]) |
                 ((time[:, None] == time[None, :]) & ~event[None, :]))
    ip, jp = np.where(comparable)
    difference = matrix[ip] - matrix[jp]
    agreement = (difference > 1e-8) + .5 * (np.abs(difference) <= 1e-8)
    result = np.full((B, len(names), 1+len(times)), np.nan)
    for b in range(B):
        index = rng.integers(n, size=n)
        counts = np.bincount(index, minlength=n)
        pairw = counts[ip] * counts[jp]
        result[b,:,0] = np.sum(pairw[:,None]*agreement, axis=0)/pairw.sum()
        for j in range(len(names)):
            for k,h in enumerate(times):
                result[b,j,k+1] = roc_auc(time[index], event[index], matrix[index,j], h)
    rows, draws = [], []
    labels = [('C-index',np.nan)] + [('AUC',h) for h in times]
    for j,name in enumerate(names):
        points = [cindex(time,event,matrix[:,j])] + [roc_auc(time,event,matrix[:,j],h) for h in times]
        for k,(metric,h) in enumerate(labels):
            v=result[:,j,k];good=v[np.isfinite(v)]
            lo,hi = np.quantile(good,[.025,.975]) if len(good) else [np.nan,np.nan]
            rows.append(dict(predictor=name,metric=metric,time_days=h,estimate=points[k],ci_lower=lo,ci_upper=hi,
                n=n,events=int(event.sum()),bootstrap_B=B,bootstrap_valid=len(good),ci_method='paired sample-row percentile bootstrap; fixed scores',
                cases_before_t=int(((time<h)&event).sum()) if np.isfinite(h) else np.nan,
                controls_after_t=int((time>h).sum()) if np.isfinite(h) else np.nan,
                fragile_horizon=bool((time>h).sum()<10) if np.isfinite(h) else False))
            draws.extend(dict(replicate=b+1,predictor=name,metric=metric,time_days=h,value=v[b]) for b in range(B))
    return pd.DataFrame(rows), pd.DataFrame(draws), result
