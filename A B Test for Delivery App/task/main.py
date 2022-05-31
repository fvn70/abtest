import pandas as pd
import scipy.stats as st
import statsmodels.stats.power as ssp


def levene_test():
    df = pd.read_csv('aa_test.csv')
    s1 = df.iloc[:, 0]
    s2 = df.iloc[:, 1]

    print("Levene's test")
    w, p = st.levene(s1, s2, center='mean')
    w = round(w, 3)
    if p > 0.05:
        pv = ">"
        h1 = 'no'
        eq = 'yes'
    else:
        pv = "<="
        h1 = 'yes'
        eq = 'no'
    print(f"W = {w}, p-value {pv} 0.05")
    print(f"Reject null hypothesis: {h1}")
    print(f"Variances are equal: {eq}")


def t_test():
    global s1, s2
    print("\nT-test")
    t, p = st.ttest_ind(s1, s2, equal_var=(eq == 'yes'))
    t = round(t, 3)
    if p > 0.05:
        pv = ">"
        h1 = 'no'
        eq = 'yes'
    else:
        pv = "<="
        h1 = 'yes'
        eq = 'no'
    print(f"t = {t}, p-value {pv} 0.05")
    print(f"Reject null hypothesis: {h1}")
    print(f"Means are equal: {eq}")

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

power_analysis()
