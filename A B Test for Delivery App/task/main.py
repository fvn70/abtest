import pandas as pd
import scipy.stats as st
import statsmodels.stats.power as ssp
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt


def hyp_by_p(sv, p, stat, text):
    if p > 0.05:
        pv = ">"
        h1 = 'no'
        eq = 'yes'
    else:
        pv = "<="
        h1 = 'yes'
        eq = 'no'
    print(f"{stat} = {sv}, p-value {pv} 0.05")
    print(f"Reject null hypothesis: {h1}")
    print(f"{text}: {eq}")
    return eq == 'yes'


def levene_test():
    df = pd.read_csv('aa_test.csv')
    s1 = df.iloc[:, 0]
    s2 = df.iloc[:, 1]
    print("Levene's test")
    w, p = st.levene(s1, s2, center='mean')
    w = round(w, 3)
    hyp_by_p(w, p, 'W', 'Variances are equal')

def t_test(eq):
    global s1, s2
    print("\nT-test")
    t, p = st.ttest_ind(s1, s2, equal_var=(eq == 'yes'))
    t = round(t, 3)
    hyp_by_p(t, p, 't', 'Means are equal')

def power_analysis():
    df = pd.read_csv('ab_test.csv')
    s1 = df.order_value[df.group == 'Control']
    s2 = df.order_value[df.group == 'Experimental']
    n1 = s1.count()
    n2 = s2.count()

    ef = 0.2
    alpha = 0.05
    power = 0.8
    ratio = n2 / n1
    n0 = ssp.tt_ind_solve_power(effect_size=ef, alpha = alpha, power = power, ratio = ratio)
    n0 = (int(n0 / 100) + 1) * 100
    print(f"Sample size: {n0}")

    n1 = min(n0, n1)
    n2 = min(n0, n1 * ratio)
    print(f"\nControl group: {n1}")
    print(f"Experimental group: {n2}")

def sessions(df):
    s_a = df.date[df.group == 'Control'].value_counts().sort_index()
    s_b = df.date[df.group == 'Experimental'].value_counts().sort_index()
    days = [d[-2:] for d in s_a.index]

    x_axis = np.arange(len(days))

    plt.figure(figsize=(10, 6))
    plt.bar(x_axis-0.2, s_a, width=0.4, label='Control')
    plt.bar(x_axis+0.2, s_b, width=0.4, label='Experimental')

    plt.xticks(x_axis, days)
    plt.xlabel(dt.strftime(pd.to_datetime(df.date[0]), '%B'), fontsize=14)
    plt.ylabel('Number of sessions', fontsize=14)
    plt.legend()

    plt.show()

def hist(df, col, xlabel):
    fig, ax = plt.subplots(1, 2, sharex=True, sharey=True)
    ax = df.hist(column=col, by='group', ax=ax)
    fig.supylabel('Frequency')
    fig.supxlabel(xlabel, fontsize=14)
    plt.show()

def stat(df):
    s = remove_outliers(df).order_value
    m = s.mean()
    sd = s.values.std(ddof=0)
    mx = s.max()
    print(f"Mean: {round(m, 2)}")
    print(f"Standard deviation: {round(sd, 2)}")
    print(f"Max: {round(mx, 2)}")

def mwu_test(df):
    d = remove_outliers(df)
    s_a = d.order_value[d.group == 'Control']
    s_b = d.order_value[d.group == 'Experimental']
    u1, p = st.mannwhitneyu(s_a, s_b)
    print("Mann-Whitney U test")
    hyp_by_p(u1, p, 'U1', 'Distributions are same')

def remove_outliers(df):
    q1 = df.order_value.quantile(0.99)
    q2 = df.session_duration.quantile(0.99)
    return df[(df.order_value < q1) & (df.session_duration < q2)]

def normalize(df):
    d0 = remove_outliers(df)
    d = pd.DataFrame({'group': d0.group, 'log_order_value': np.log(d0.order_value)})
    d.plot(y='log_order_value', kind='hist', bins=20)
    plt.ylabel('Frequency')
    plt.xlabel('Log order value', fontsize=14)
    plt.show()

    print("Levene's test")
    s_a = d.log_order_value[d.group == 'Control']
    s_b = d.log_order_value[d.group == 'Experimental']

    w, p = st.levene(s_a, s_b)
    w = round(w, 3)
    eq = hyp_by_p(w, p, 'W', 'Variances are equal')
    print("\nT-test")
    t, p = st.ttest_ind(s_a, s_b, equal_var=eq)
    t = round(t, 3)
    hyp_by_p(t, p, 't', 'Means are equal')



df = pd.read_csv('ab_test.csv')

# stage2
# power_analysis()

# stage3
# sessions(df)
# hist(df, 'order_value', 'Order value')
# hist(df, 'session_duration', 'Session duration')
# stat(df)

# stage4
# mwu_test(df)

normalize(df)
