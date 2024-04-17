using DelimitedFiles

cotes::Array{Int64,2} = readdlm("cotes.dat")

function main()
    function f(x)
        return 1/(x+1)
    end

    for i in UnitRange(1,8)
        display(newtoncotes(f,0,1,i))
    end

    display("Trapecio Compuesto:")
    display(trapeciocompuesto(f,0,1,100000))
    
    display("Simpson Compuesto:")
    display(simpsoncompuesto(f,0,1,100000))

    display("Romberg")
    display(romberg(f,0,1,30))

    display("Derivación")
    display(derivadadefinicion(f,1,0.0001))
    display(derivadataylor(f,1,0.0001))
    display(derivadarichardson(f,1,0.0001))
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
    if n_intervalos%2 != 0
        return "error, n ha de ser par"
    end
    
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

function romberg(func::Function, a, b, p::Int64)
    if p == 0
        return ((b-a)/2)*(func(a)+func(b))
    else
        suma = 0
        h2p = (b-a)/(2^p)

        for i in range(1,2^(p-1))
            suma += func(a+(2*i-1)*h2p)
        end

        return 0.5*romberg(func,a,b,p-1) + h2p*suma
    end
end

function derivadadefinicion(func::Function, x, h)
    return (func(x+h)-func(x))/h
end

function derivadataylor(func::Function, x, h)
    return (func(x+h)-func(x-h))/(2*h)
end

function derivadarichardson(f::Function, x, h)
    return (4.0/3.0)*((f(x+h/2)-f(x-h/2))/h) - (1.0/3.0)*((f(x+h)-f(x-h))/(2*h))
end

main()
