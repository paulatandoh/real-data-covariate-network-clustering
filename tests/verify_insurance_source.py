"""Optional source-data check; requires pandas. Run from the project root."""
from pathlib import Path
import csv
import pandas as pd
root=Path(__file__).resolve().parents[1]
df=pd.read_stata(root/'data/raw/0422allinforawnet.dta',convert_categoricals=False)
df=df.loc[(df.village.str.lower()=='dukou')&(df.network_village.str.lower()=='dukou')]
df=df.loc[df.id.notna()&df.network_id.notna()&(df.network_id!=99)&(df.id!=df.network_id)]
t=df.loc[df.takeup_survey.isin([0,1]),['network_id','takeup_survey']].drop_duplicates()
assert t.groupby('network_id').takeup_survey.nunique().max()==1
eligible=set(t.network_id)
df=df.loc[df.id.isin(eligible)&df.network_id.isin(eligible)]
edges={tuple(sorted((int(a),int(b))))for a,b in zip(df.id,df.network_id)}
cov=list(csv.DictReader(open(root/'data/insurance/insurance_covariate.csv')))
ids=[int(r['original_farmer_id'])for r in cov]
lookup=dict(zip(t.network_id,t.takeup_survey))
assert all(int(r['insured'])==lookup[int(r['original_farmer_id'])]for r in cov)
with open(root/'data/insurance/insurance_adjacency.csv')as f:rows=list(csv.reader(f))[1:]
saved={tuple(sorted((ids[i],ids[j])))for i,row in enumerate(rows)for j,v in enumerate(row[1:])if int(v)==1 and i<j}
assert len(edges)==667 and len({v for e in edges for v in e})==250 and saved==edges
print('PASS: original Stata reconstruction equals the supplied insurance graph and covariate lookup (250 nodes / 667 edges).')
