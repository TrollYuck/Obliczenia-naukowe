#Tomasz Niedziałek 279754
using LinearAlgebra, Printf

function matcond(n::Int, c::Float64)
# Function generates a random square matrix A of size n with
# a given condition number c.
# Inputs:
#	n: size of matrix A, n>1
#	c: condition of matrix A, c>= 1.0
#
# Usage: matcond(10, 100.0)
#
# Pawel Zielinski
        if n < 2
         error("size n should be > 1")
        end
        if c< 1.0
         error("condition number  c of a matrix  should be >= 1.0")
        end
        (U,S,V)=svd(rand(n,n))
        return U*diagm(0 =>[LinRange(1.0,c,n);])*V'
end

function hilb(n::Int)
# Function generates the Hilbert matrix  A of size n,
#  A (i, j) = 1 / (i + j - 1)
# Inputs:
#	n: size of matrix A, n>=1
#
#
# Usage: hilb(10)
#
# Pawel Zielinski
        if n < 1
         error("size n should be >= 1")
        end
        return [1 / (i + j - 1) for i in 1:n, j in 1:n]
end

function runhilbert(n::Int)
    # @assert n >= 1
    println("----Hilbert----")
    for i in 1:n
        A = hilb(i)
        x = ones(Float64, i)
        b = A * x

        xGauss = A \ b
        xInv = inv(A) * b

        xErrorGauss = norm(xGauss - x) / norm(x)
        xErrorInv = norm(xInv - x) / norm(x)

        println("__________________________________")
        @printf("Size: %d x %d | Condition: %e | Rank: %d \n", i, i, cond(A), rank(A));
        @printf("Relative error for Gauss: %e | Relative error for inv: %e\n", xErrorGauss, xErrorInv)
    end
    println("\n\n")
end

function runrandom(n::Array{Int},c::Array{Float64})
    println("----Random----")
    for i in n
        for j in c
            A = matcond(i,j)
            x = ones(Float64, i)
            b = A * x 

            xGauss = A \ b 
            xInv = inv(A) * b 

            xErrorGauss = norm(xGauss - x) / norm(x)
            xErrorInv = norm(xInv - x) / norm(x)

            println("__________________________________") 
            @printf("Size: %d x %d | Condition: %e | Rank: %d \n", i, i, j, rank(A));
            @printf("Relative error for Gauss: %e | Relative error for inv: %e\n", xErrorGauss, xErrorInv)

        end
    end
    println("\n\n")
end

sizes = [5,10,20]
cs = [Float64(1), Float64(10), Float64(10^3), Float64(10^7), Float64(10^12), Float64(10^16)]

runhilbert(20)
runrandom(sizes, cs)