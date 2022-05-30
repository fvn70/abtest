import pandas as pd
import scipy.stats as st

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
