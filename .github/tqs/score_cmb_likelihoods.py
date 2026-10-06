#!/usr/bin/env python3
import argparse, json, os, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

TCMB_K = 2.7255
UK2 = (TCMB_K*1e6)**2

def read_class_dl(path):
    arr=np.loadtxt(path)
    ell=arr[:,0].astype(int)
    lmax=int(ell.max())
    out={k:np.zeros(lmax+1) for k in ["tt","ee","te","bb"]}
    # CLASS cl_lensed.dat stores dimensionless D_l=l(l+1)C_l/2pi
    out["tt"][ell]=arr[:,1]*UK2
    out["ee"][ell]=arr[:,2]*UK2
    out["te"][ell]=arr[:,3]*UK2
    if arr.shape[1]>4:
        out["bb"][ell]=arr[:,4]*UK2
    return out,lmax

def planck_score(class_file, planck_repo):
    sys.path.insert(0,str(Path(planck_repo).resolve()))
    from planck_lite_py import PlanckLitePy
    cl,lmax=read_class_dl(class_file)
    if lmax < 2508:
        raise RuntimeError(f"Planck 2018 Plik-lite needs ell>=2508, got {lmax}")
    like=PlanckLitePy(data_directory=str(Path(planck_repo)/"data"),
                      year=2018,spectra="TTTEEE",use_low_ell_bins=False)
    logl=float(like.loglike(cl["tt"],cl["te"],cl["ee"],ellmin=0))
    return {"loglike":logl,"chi2_relative_normalization":-2*logl,
            "ell_min":30,"ell_max_TT":2508,"ell_max_TE_EE":1996}

def locate_act_fits(root):
    hits=list(Path(root).rglob("dr6_data_cmbonly.fits"))
    if not hits:
        raise FileNotFoundError(f"Could not find dr6_data_cmbonly.fits below {root}")
    return hits[0]

def act_prepare(data_file):
    import sacc
    inp=sacc.Sacc.load_fits(str(data_file))
    cuts={"TT":(600,6500),"TE":(600,6500),"EE":(600,6500)}
    meta=[]; cull=[]; idxmax=0
    pol_dt={"t":"0","e":"e","b":"b"}
    for pol in ["TT","TE","EE"]:
        p1,p2=pol.lower()
        dt=f"cl_{pol_dt[p1]}{pol_dt[p2]}"
        for tr1,tr2 in inp.get_tracer_combinations(dt):
            ls,mu,ind=inp.get_ell_cl(dt,tr1,tr2,return_ind=True)
            mask=(ls>=cuts[pol][0])&(ls<=cuts[pol][1])
            if not np.all(mask):
                cull.append(ind[~mask])
            if np.any(mask):
                window=inp.get_bandpower_windows(ind[mask])
                meta.append(dict(pol=pol.lower(),idx=ind[mask],spec=mu[mask],
                                 win=window.weight.T,ells=window.values.astype(int)))
                idxmax=max(idxmax,int(np.max(ind)))
    data=np.zeros(idxmax+1)
    for m in meta:
        data[m["idx"]]=m["spec"]
    cov=inp.covariance.covmat.copy()
    for cc in cull:
        cov[cc,:]=0.; cov[:,cc]=0.; cov[cc,cc]=1e10
    return data,np.linalg.inv(cov),meta

def act_score(class_file, act_data_root):
    cl,lmax=read_class_dl(class_file)
    if lmax < 9000:
        raise RuntimeError(f"ACT DR6 CMB-only requests theory to ell=9000, got {lmax}")
    data,icov,meta=act_prepare(locate_act_fits(act_data_root))
    def chi(x):
        Aact,Pact=x
        pred=np.zeros_like(data)
        for m in meta:
            pol=m["pol"]; ells=m["ells"]
            if np.max(ells)>lmax:
                return 1e100
            dat=cl[pol][ells]/(Aact*Aact)
            if pol[0]=="e": dat=dat/Pact
            if pol[1]=="e": dat=dat/Pact
            pred[m["idx"]]=m["win"]@dat
        d=data-pred
        return float(d@icov@d)
    opt=minimize(chi,[1.,1.],method="Nelder-Mead",
                 options={"maxiter":2000,"xatol":1e-8,"fatol":1e-5})
    # enforce official flat bounds with a bounded polish
    opt2=minimize(chi,np.clip(opt.x,[.5,.9],[1.5,1.1]),method="L-BFGS-B",
                  bounds=[(.5,1.5),(.9,1.1)],options={"maxiter":2000,"ftol":1e-12})
    return {"chi2":float(opt2.fun),"loglike":float(-0.5*opt2.fun),
            "A_act":float(opt2.x[0]),"P_act":float(opt2.x[1]),
            "ell_min":600,"ell_max":6500,"theory_lmax_required":9000}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--lcdm",required=True)
    ap.add_argument("--tqs")
    ap.add_argument("--planck-repo",required=True)
    ap.add_argument("--act-data",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    res={"lcdm":{},"tqs":{},"delta_chi2":{}}
    for name,path in [("lcdm",a.lcdm),("tqs",a.tqs)]:
        if not path or not Path(path).exists():
            res[name]={"status":"spectra_unavailable"}
            continue
        res[name]["status"]="ok"
        res[name]["planck2018_pliklite"]=planck_score(path,a.planck_repo)
        res[name]["act_dr6_cmbonly"]=act_score(path,a.act_data)
    if res["lcdm"].get("status")=="ok" and res["tqs"].get("status")=="ok":
        for key in ["planck2018_pliklite","act_dr6_cmbonly"]:
            res["delta_chi2"][key]=(
                -2*res["tqs"][key]["loglike"]+2*res["lcdm"][key]["loglike"]
            )
        res["delta_chi2"]["combined"]=sum(res["delta_chi2"].values())
    Path(a.output).write_text(json.dumps(res,indent=2))
    print(json.dumps(res,indent=2))
if __name__=="__main__":
    main()
