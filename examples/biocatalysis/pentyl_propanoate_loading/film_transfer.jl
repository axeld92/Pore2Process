module LoadingFilmTransfer

import ..LorenzoPackedBed: rate, rate_terms
import ..BeadFilm

export LoadingAdjustedFilm, inlet_balance

"""
Apply an independently derived mass-specific rate multiplier inside
the existing liquid-film balance. The axial solver still receives the
measured catalyst concentration; film area and all geometric inventory
quantities are unchanged. An enhancement does not imply eta > 1.
"""
struct LoadingAdjustedFilm
    film::BeadFilm.FilmKinetics
    relative_rate::Float64
    function LoadingAdjustedFilm(film, relative_rate)
        isfinite(relative_rate) && relative_rate > 0 ||
            throw(ArgumentError("Use a finite positive rate ratio."))
        new(film, relative_rate)
    end
end

rate_terms(model::LoadingAdjustedFilm, args...) = rate_terms(model.film, args...)

function rate(model::LoadingAdjustedFilm, terms, acid0, alcohol0, catalyst_g_L, conversion)
    # The published local rate is linear in its catalyst argument. The
    # factor acts here solely on the chemical source term before solving
    # for surface concentration, not on the catalyst inventory.
    rate(model.film, terms, acid0, alcohol0,
         catalyst_g_L*model.relative_rate, conversion)
end

inlet_balance(model::LoadingAdjustedFilm, T, catalyst_g_L; acid0=3.3, alcohol0=6.6) =
    BeadFilm.inlet_balance(model.film, T, catalyst_g_L*model.relative_rate; acid0, alcohol0)

end
