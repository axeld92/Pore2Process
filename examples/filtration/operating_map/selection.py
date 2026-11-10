"""Select only among actually solved grid points under declared policies."""
from itertools import product

def selections(cases,p):
    rows=[]
    supports=sorted({c['request']['support'] for c in cases})
    for support,flux,fixed,mapping,fraction,factor,policy in product(supports,p['process']['fluxes_LMH'],
            p['process']['fixed_multipliers'],p['capacity']['mappings'],p['capacity']['fractions'],
            p['capacity']['transfer_factors'],p['capacity']['policies']):
        matches=[];evaluated=0
        for c in cases:
            q=c['request']
            if q.get('viscosity_multiplier',1.)!=1.:continue
            if (q['support'],q['flux_LMH'],q['fixed_multiplier'])!=(support,flux,fixed):continue
            d=next(v for v in c['capacity'] if (v['mapping'],v['capacity_fraction'],v['transfer_factor'])==(mapping,fraction,factor))
            if d['status']!='evaluated':continue
            evaluated+=1
            permitted=d['nominal_admissible'] if policy=='linear_nominal' else d['linear_reading_admissible'] if policy=='linear_lower_reading' else d['monotonic_classification']=='admissible'
            if permitted:matches.append(c)
        best=min(matches,key=lambda c:c['validation']['specific_shaft_energy_kWh_m3_reference_volume'],default=None)
        rows.append({'support':support,'mean_flux_LMH':flux,'fixed_multiplier':fixed,'mapping':mapping,
            'capacity_fraction':fraction,'transfer_factor':factor,'policy':policy,'candidates':len(matches),
            'evaluated_source_mappings':evaluated,'selected_case_id':None if best is None else best['id'],
            'selected_Re':None if best is None else best['request']['inlet_Re'],
            'specific_shaft_energy_kWh_m3':None if best is None else best['validation']['specific_shaft_energy_kWh_m3_reference_volume'],
            'claim':'Minimum on evaluated discrete grid under stated transfer and capacity assumptions; not a validated physical optimum.'})
    return rows
