"""Usage: python scripts/tables/prediction_table.py CONFIG.json NEW_OUTPUT"""
import csv,json,math,hashlib,sys
from pathlib import Path
from statistics import mean
if len(sys.argv)!=3:raise SystemExit(__doc__)
config=Path(sys.argv[1]).resolve();cfg=json.loads(config.read_text());out=Path(sys.argv[2]);paths={}
for panel,key in [('A','nested_summary'),('B','loeo_metrics')]:
 entry=cfg['inputs'][key];path=(config.parent/entry['path']).resolve()
 if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Input hash mismatch')
 paths[panel]=path
panels=[]
headers=['Environment','N','Median_Pearson_r','All_Pearson_r','Median_RMSE_t_ha','All_RMSE_t_ha','RMSE_Reduction_Percent']
for panel,p in paths.items():
    source=list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
    source=[x for x in source if x['Scenario'] in ('Median','all') and (panel=='A' or x['Scaling']=='WITHIN')]
    assert len(source)==22
    key='Env' if panel=='A' else 'Heldout'
    ncol='N_Genotypes' if panel=='A' else 'n_test'
    rcol='Mean_Pearson' if panel=='A' else 'Pearson'
    ecol='Mean_RMSE' if panel=='A' else 'RMSE'
    keys=[(x[key],x['Scenario']) for x in source]
    assert len(set(keys))==22
    envs=sorted(set(x[key] for x in source));assert len(envs)==11
    lookup={k:x for k,x in zip(keys,source)}
    rows=[]
    for env in envs:
        m,a=[lookup[(env,s)] for s in ['Median','all']]
        assert int(m[ncol])==int(a[ncol])
        if panel=='A': assert int(m['N_Repeats'])==int(a['N_Repeats'])==5
        else: assert m['SourceStatus']==a['SourceStatus']=='historical_not_refit'
        nums=[float(m[rcol]),float(a[rcol]),float(m[ecol]),float(a[ecol])]
        assert all(math.isfinite(x) for x in nums) and nums[2]>0 and nums[3]>0
        rows.append([env,int(m[ncol]),*nums,100*(nums[2]-nums[3])/nums[2]])
    overall=['Equal-weight mean',None]+[mean(r[j] for r in rows) for j in range(2,7)]
    file='Table_S3A_Within_Environment.csv' if panel=='A' else 'Table_S3B_LOEO_Within_Scaling.csv'
    panels.append(dict(panel=panel,file=file,rows=rows,overall=overall,total_environment_genotype_pairs=sum(r[1] for r in rows),source=str(p),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source_status='V15 nested primary' if panel=='A' else '22 historical_not_refit WITHIN Median/all models',count_both_improved=sum(r[3]>r[2] and r[5]<r[4] for r in rows)))

assert [r[0:2] for r in panels[0]['rows']]==[r[0:2] for r in panels[1]['rows']]
assert abs(panels[0]['overall'][2]-.495498164)<1e-8
assert abs(panels[0]['overall'][3]-.571019194)<1e-8
assert abs(panels[0]['overall'][6]-4.952196)<1e-5
assert panels[0]['count_both_improved']==10
assert panels[0]['total_environment_genotype_pairs']==6398
out.mkdir(exist_ok=False,parents=True)
for p in panels:
    dest=out/p['file'];assert not dest.exists(),dest
    with dest.open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f);writer.writerow(headers);writer.writerows(p['rows']);writer.writerow(p['overall'])
