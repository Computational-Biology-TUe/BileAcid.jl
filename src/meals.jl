# Meals and the postprandial (gastro-intestinal) reflex.

"""
    gi_reflex(β, δ, t)

Postprandial reflex at time `t` after the start of a meal: a transient multiplier of a
basal rate that rises from 1 at `t = 0` to its peak `1 + δ` at `t = β` and returns to 1
afterwards. It is 1 before the meal.
"""
@inline function gi_reflex(β, δ, t)
    r = 1 + δ * (t / β) * exp(1 / 2 - t^2 / (2 * β^2))
    return max(r, one(r))
end

"""
    postprandial_reflex(m::BileAcidModel, t) -> (gallbladder, gut)

Multipliers of gallbladder emptying and intestinal transit at time `t`. As in the MATLAB
model only the most recent meal counts; before the first meal both are 1.
"""
@inline function postprandial_reflex(m::BileAcidModel, t)
    i = searchsortedlast(m.meal_times, t)
    since_meal = i == 0 ? zero(t) : t - m.meal_times[i]
    gb = gi_reflex(m.gb_reflex..., since_meal)
    gut = gi_reflex(m.gut_reflex..., since_meal)
    return i == 0 ? (one(gb), one(gut)) : (gb, gut)
end

"""
    MealSchedule(; fasting_days = 20, meal_days = 35, gaps = (6, 6, 12))
    MealSchedule(meal_times, tspan)

When meals start and how long to simulate (both in minutes).

The keyword form is the protocol of the MATLAB implementation (`simulate_human_day`):
`fasting_days` days without meals, to build up the bile acid pool, followed by
`meal_days` days with meals separated by `gaps` hours, repeated daily. As in the MATLAB
implementation, the fasting period starts with a postprandial response at `t = 0`.
"""
struct MealSchedule
    meal_times::Vector{Float64}
    tspan::Tuple{Float64, Float64}
end

function MealSchedule(; fasting_days::Real = 20, meal_days::Integer = 35, gaps = (6, 6, 12))
    fasting_end = 24 * 60 * float(fasting_days)
    intervals = repeat(60 .* collect(Float64, gaps), meal_days)
    meal_starts = fasting_end .+ cumsum([0.0; intervals])[1:(end - 1)]
    meal_times = unique!([0.0; meal_starts])        # fasting_days = 0: first meal at t = 0
    return MealSchedule(meal_times, (0.0, fasting_end + sum(intervals)))
end

"""
    last_day(s::MealSchedule)

Time span of the final 24 hours of the schedule.
"""
last_day(s::MealSchedule) = (s.tspan[2] - 24 * 60, s.tspan[2])
