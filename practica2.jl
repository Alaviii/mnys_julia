using DelimitedFiles

function f(x)
    return 1/(x+1)
end

function newtoncotes(func::Function, lim_inf, lim_sup, n_newton_cotes::Int64)
    integral::Float64 = 0

    h::Float64 = (lim_sup-lim_inf)/n_newton_cotes

    for i in UnitRange(1,n_newton_cotes+1)
        integral += (cotes[n_newton_cotes,i]/cotes[n_newton_cotes,10])*func(lim_inf+(i-1)*h)
    end
    
    return (lim_sup-lim_inf)*integral
end

function trapeciocompuesto(func::Function, lim_inf, lim_sup, n_intervalos::Int64)
    h::Float64 = (lim_sup-lim_inf)/n_intervalos

    integral::Float64 = (h/2)*(func(lim_inf)+func(lim_sup))

    suma::Float64 = 0

    for i in UnitRange(1,n_intervalos-1)
        suma += func(lim_inf + i*h)
    end
    
    return integral + h*suma
end

function simpsoncompuesto(func::Function, lim_inf, lim_sup, n_intervalos::Int64)
    h::Float64 = (lim_sup-lim_inf)/n_intervalos

    suma_pares::Float64 = 0
    suma_impares::Float64 = 0

    i::Int64 = 1

    while i<n_intervalos
        suma_impares += func(lim_inf + i*h)
        i+= 2
    end

    i = 2

    while i<n_intervalos
        suma_pares += func(lim_inf + i*h)
        i+= 2
    end
    
    return (h/3)*(func(lim_inf)+func(lim_sup)+4*suma_impares+2*suma_pares)
end


for i in UnitRange(1,8)
    display(newtoncotes(f,0,1,i))
end

cotes::Array{Int64,2} = readdlm("cotes.dat")

display(trapeciocompuesto(f,0,1,10000))

display(simpsoncompuesto(f,0,1,100000))