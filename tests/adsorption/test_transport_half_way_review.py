"""Independent coverage oracle catches absent/spurious reflection links."""
from itertools import product
import numpy as np
import pytest
from examples.adsorption.transport_root_cause.flux_audit.half_way.independent_review import descriptor,expected_links,coverage,mask


def stock_contract():
    vectors=[(0,0,0)]+[x for x in product((-1,0,1),repeat=3) if x!=(0,0,0)]
    rows=[]
    for p,v in enumerate(vectors):
        opposite=vectors.index(tuple(-x for x in v));weight=np.prod([2/3 if x==0 else 1/6 for x in v])
        rows.append((p,opposite,*v,weight))
    return np.array(rows,dtype=[('p',int),('opposite',int),('cx',int),('cy',int),('cz',int),('weight',float)])


def channel():
    materials=np.ones((4,6,4),dtype=int);materials[:,0,:]=2;materials[:,-1,:]=2
    directions=descriptor(stock_contract())[0];expected,_=expected_links(materials,directions,(True,False,True))
    links=np.array([(*key,.5) for key in sorted(expected)],dtype=[('i',int),('j',int),('k',int),('p',int),('q',float)])
    return materials,directions,links


def test_flat_links_include_all_nine_wall_directions_and_periodic_seams():
    materials,directions,links=channel();report=coverage(materials,directions,links,(True,False,True))
    assert report['expected_fluid_solid_links']==2*4*4*9
    assert report['expected_periodic_seam_links']>0
    assert report['covered_periodic_seam_links']==report['expected_periodic_seam_links']
    assert report['passed']


def test_missing_periodic_link_is_detected():
    materials,directions,links=channel();_,seams=expected_links(materials,directions,(True,False,True))
    bad=next(iter(seams));keep=[tuple(int(row[k]) for k in ('i','j','k','p'))!=bad for row in links]
    report=coverage(materials,directions,links[np.array(keep)],(True,False,True))
    assert not report['passed'] and report['missing_links']==1
    assert report['covered_periodic_seam_links']==report['expected_periodic_seam_links']-1


def test_fluid_fluid_reflection_and_non_half_distance_are_rejected():
    materials,directions,links=channel();p=list(map(tuple,directions)).index((1,0,0))
    extra=np.array([(1,3,1,p,.5)],dtype=links.dtype)
    assert coverage(materials,directions,np.concatenate([links,extra]),(True,False,True))['unexpected_links']==1
    links=links.copy();links['q'][0]=0.
    report=coverage(materials,directions,links,(True,False,True))
    assert not report['passed'] and report['non_halfway_links']==1


def test_descriptor_opposite_and_incomplete_mask_fail():
    table=stock_contract();table['opposite'][1]=1
    with pytest.raises(ValueError):descriptor(table)
    incomplete=np.array([(0,0,0,1),(1,1,1,2)],dtype=[('i',int),('j',int),('k',int),('material',int)])
    with pytest.raises(ValueError):mask(incomplete)
