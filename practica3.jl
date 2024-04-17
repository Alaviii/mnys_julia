using DelimitedFiles

function main()
    matrix::Array{Float64, 2} = readdlm("matrizin.dat")

    writedlm("matrizout.dat",matrix)
end

main()