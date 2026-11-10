"""Engineering screens, separate from solver and software success."""
import numpy as np

from .data import evidence


def batch_metrics(result,data):
    t=data['time_s'];q=data['loading_kg_kg']
    pred=np.interp(t,result['time_s'],result['mean_loading_kg_kg'])
    e=evidence()['figure_2'];dt=e['reading_sensitivity_time_s'];dq=e['reading_sensitivity_loading_kg_kg']
    corners=[]
    for sx in [-1,1]:
        for sy in [-1,1]:
            p=np.interp(np.maximum(t+sx*dt,0),result['time_s'],result['mean_loading_kg_kg'])
            corners.append(float(np.sqrt(np.mean((p-(q+sy*dq))**2)))*1000)
    rmse=float(np.sqrt(np.mean((pred-q)**2)))*1000
    return dict(RMSE_mg_g=rmse,maximum_absolute_error_mg_g=float(np.max(np.abs(pred-q)))*1000,
        reading_sensitivity_RMSE_mg_g=[min(corners),max(corners)],RMSE_screen_passed=rmse<=3.,
        role=data['role'],diameter_m=data['diameter_m'],points=len(t),
        normalized_mass_residual=result['normalized_mass_residual'],
        normalized_interface_flux_mismatch=result['normalized_interface_flux_mismatch'])


def column_passed(m):
    return m['RMSE_screen_passed'] and m['t50_screen_passed'] and m['t05_screen_passed'] is not False


def physical_status(endpoint,batch_records,column_records):
    batch_ok=endpoint['source_endpoint_screen_passed'] and all(r['metrics']['RMSE_screen_passed'] for r in batch_records)
    by_id={r['id']:r for r in column_records}
    primary=all(column_passed(by_id[f'primary_{s}']['metrics']) for s in ['calibration','holdout'])
    calibration=by_id['anchored_calibration']['metrics']
    anchored=calibration['RMSE_screen_passed'] and calibration['t50_screen_passed'] and column_passed(by_id['anchored_holdout']['metrics'])
    return dict(source_equilibrium_audit_passed=endpoint['source_endpoint_screen_passed'],
        batch_particle_size_comparison_passed=batch_ok,unfitted_column_comparison_passed=primary,
        one_anchor_column_comparison_passed=anchored,
        full_transport_comparison_passed=batch_ok and (primary or anchored),
        experimental_pressure_validation=False)
