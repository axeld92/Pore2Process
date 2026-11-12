"""Peak localization must retain diagonal/periodic wall proximity and fluid IDs."""
import csv
import json

import pytest

from examples.adsorption.transport_root_cause.study import peak_analysis


def test_peak_localization_uses_fluid_extrema_and_periodic_diagonal_links(tmp_path):
    # A wall at (0,0,0), a diagonal fluid maximum at the periodic corner,
    # and an interior minimum distinguish link adjacency from face-only tests.
    n=5
    with (tmp_path/'peak.csv').open('w',newline='') as stream:
        out=csv.writer(stream)
        out.writerow(['i','j','k','material','x_m','y_m','z_m','concentration_kg_m3','ux_m_s','uy_m_s','uz_m_s'])
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    c=.5
                    if (i,j,k)==(0,0,0):c=9. # numerical wall is not a fluid extremum
                    if (i,j,k)==(4,4,4):c=.9
                    if (i,j,k)==(2,2,2):c=.1
                    out.writerow([i,j,k,2 if (i,j,k)==(0,0,0) else 1,i+.5,j+.5,k+.5,c,0,0,0])
    (tmp_path/'peak.json').write_text(json.dumps(dict(step=10,time_s=.01,maximum_excursion_kg_m3=.15,initial_amplitude_kg_m3=.25)))
    r=peak_analysis(tmp_path)
    assert r['maximum']['ijk']==[4,4,4]
    assert r['maximum']['touches_solid_link']
    assert not r['minimum']['touches_solid_link']
    assert r['fraction_fluid_outside_engineering_bounds']==pytest.approx(2/124)
    assert r['fraction_violators_touching_solid_link']==.5
    assert r['fraction_fluid_touching_solid_link']==pytest.approx(26/124)
