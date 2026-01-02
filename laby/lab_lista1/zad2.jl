#Tomasz Niedziałek 279754

println("Kahan half eps = ", Float16(Float16(3)*(Float16(4)/Float16(3)-Float16(1))-Float16(1)), ", Julia = ", eps(Float16))
println("Kahan single eps = ", Float32(Float32(3)*(Float32(4)/Float32(3)-Float32(1))-Float32(1)), ", Julia = ", eps(Float32))
println("Kahan double eps = ", Float64(Float64(3)*(Float64(4)/Float64(3)-Float64(1))-Float64(1)), ", Julia = ", eps(Float64))