using DelimitedFiles

cotes::Array{Int64,2} = readdlm("cotes.dat")

writedlm("test.dat",cotes)

function f(x::Float64)
    return 1/(x+1)
end

function newtoncotessimple(func::Function, lim_inf, lim_sup, n_newton_cotes::Int64)
    integral::Float64 = 0

    h::Float64 = (lim_sup-lim_inf)/n_newton_cotes

    for i in UnitRange(1,n_newton_cotes+1)
        integral += (cotes[n_newton_cotes,i]/cotes[n_newton_cotes,10])*func(lim_inf+(i-1)*h)
    end
    
    return (lim_sup-lim_inf)*integral
end

for i in UnitRange(1,8)
    display(newtoncotessimple(f,0.0,1.0,i))
end