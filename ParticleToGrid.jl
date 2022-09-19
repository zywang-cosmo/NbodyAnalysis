using FFTW, Dierckx, DelimitedFiles, HDF5, Statistics


#Nearest grid point, with equal particle mass.
function NGP(posx_p, posy_p, posz_p, L, Ncell)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #end
    NGP_array = zeros(Float64, Ncell, Ncell, Ncell)
    @inbounds for i in 1:Nparticle
        nx = round(Int, posx_p[i]/H, RoundNearestTiesUp)%Ncell + 1
        ny = round(Int, posy_p[i]/H, RoundNearestTiesUp)%Ncell + 1
        nz = round(Int, posz_p[i]/H, RoundNearestTiesUp)%Ncell + 1
        NGP_array[nx, ny, nz] += 1
    end
    return NGP_array
end

#Nearest grid point, with different particle mass.
function NGP_w(posx_p, posy_p, posz_p, L, Ncell, mass_p)
    dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    if dim_p > 1
        posx_p = reshape(posx_p, Nparticle); 
        posy_p = reshape(posy_p, Nparticle); 
        posz_p = reshape(posz_p, Nparticle);
        mass_p = reshape(mass_p, Nparticle);
    end
    NGP_array = zeros(Float64, Ncell, Ncell, Ncell)
    @inbounds for i in 1:Nparticle
        mass = mass_p[i]
        nx = round(Int, posx_p[i]/H, RoundNearestTiesUp)%Ncell + 1
        ny = round(Int, posy_p[i]/H, RoundNearestTiesUp)%Ncell + 1
        nz = round(Int, posz_p[i]/H, RoundNearestTiesUp)%Ncell + 1
        NGP_array[nx, ny, nz] += 1 * mass
    end
    return NGP_array
end
# -----------------------------------------------------------------------------

#Cloud in cell, with equal particle mass.
#Interpolate the particle mass onto adjacent two grids along each dim.
#The grid points' coordinates are {0, 1, 2, 3..., Ncell-1} along each dim.
function CIC(posx_p, posy_p, posz_p, L, Ncell)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #end
    CIC_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function CIC_s(s)
        if abs(s) < 1 return 1-abs(s); else return 0
    end;end
    
    #Find the nearest two grid points around the particle along each dimension, with abs(s) < 1
    function FindxGwG_CIC!(posx::Float64, posy::Float64, posz::Float64, 
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1;
        s1 = n1 - posx;  s2 = n2 - posx;
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1;
        wx[1] = CIC_s(s1);  wx[2] = CIC_s(s2); 
        
        n1 = floor(Int, posy); n2 = n1 + 1;
        s1 = n1 - posy;  s2 = n2 - posy;
        ny[1] = (n1 + Ncell)%Ncell+1; ny[2] = (n2 + Ncell)%Ncell+1;
        wy[1] = CIC_s(s1);  wy[2] = CIC_s(s2);
        
        n1 = floor(Int, posz); n2 = n1 + 1;
        s1 = n1 - posz;  s2 = n2 - posz;
        nz[1] = (n1 + Ncell)%Ncell+1; nz[2] = (n2 + Ncell)%Ncell+1;
        wz[1] = CIC_s(s1);  wz[2] = CIC_s(s2);
    end
    
    nx, ny, nz = zeros(Int64, 2), zeros(Int64, 2), zeros(Int64, 2)
    wx, wy, wz = zeros(Float64, 2), zeros(Float64, 2), zeros(Float64, 2)
    @inbounds for i in 1:Nparticle
        FindxGwG_CIC!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx, wy, wz)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        for j1 in 1:2
            for j2 in 1:2
                for j3 in 1:2
                    CIC_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * wy[j2] * wz[j3]
    end;end;end;end
    return CIC_array
end

function CIC_w(posx_p, posy_p, posz_p, L, Ncell, mass_p)
    dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    if dim_p > 1
        posx_p = reshape(posx_p, Nparticle); 
        posy_p = reshape(posy_p, Nparticle); 
        posz_p = reshape(posz_p, Nparticle);
        mass_p = reshape(mass_p, Nparticle);
    end
    CIC_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function CIC_s(s)
        if abs(s) < 1 return 1-abs(s); else return 0
    end;end
    
    #Find the nearest two grid points around the particle along each dimension, with abs(s) < 1
    function FindxGwG_CIC!(posx::Float64, posy::Float64, posz::Float64, 
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1;
        s1 = n1 - posx;  s2 = n2 - posx;
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1;
        wx[1] = CIC_s(s1);  wx[2] = CIC_s(s2); 
        
        n1 = floor(Int, posy); n2 = n1 + 1;
        s1 = n1 - posy;  s2 = n2 - posy;
        ny[1] = (n1 + Ncell)%Ncell+1; ny[2] = (n2 + Ncell)%Ncell+1;
        wy[1] = CIC_s(s1);  wy[2] = CIC_s(s2);
        
        n1 = floor(Int, posz); n2 = n1 + 1;
        s1 = n1 - posz;  s2 = n2 - posz;
        nz[1] = (n1 + Ncell)%Ncell+1; nz[2] = (n2 + Ncell)%Ncell+1;
        wz[1] = CIC_s(s1);  wz[2] = CIC_s(s2);
    end
    
    nx, ny, nz = zeros(Int64, 2), zeros(Int64, 2), zeros(Int64, 2)
    wx, wy, wz = zeros(Float64, 2), zeros(Float64, 2), zeros(Float64, 2)
    @inbounds for i in 1:Nparticle
        mass = mass_p[i]
        FindxGwG_CIC!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx, wy, wz)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        for j1 in 1:2
            for j2 in 1:2
                for j3 in 1:2
                    CIC_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * wy[j2] * wz[j3] * mass
    end;end;end;end
    return CIC_array
end

function CIC1_w(posx_p, posy_p, posz_p, L, Ncell, mass_p)
    dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    if dim_p > 1
        posx_p = reshape(posx_p, Nparticle); 
        posy_p = reshape(posy_p, Nparticle); 
        posz_p = reshape(posz_p, Nparticle);
        mass_p = reshape(mass_p, Nparticle);
    end
    CIC_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function CIC_s(s)
        if abs(s) < 1 return 1-abs(s); else return 0
    end;end
    
    #Find the nearest two grid points around the particle along each dimension, with abs(s) < 1
    function FindxGwG_CIC!(posx::Float64, posy::Float64, posz::Float64, 
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1;
        s1 = n1 - posx;  s2 = n2 - posx;
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1;
        wx[1] = CIC_s(s1);  wx[2] = CIC_s(s2); 
        
        ny[1] = round(Int64, posy);  nz[1] = round(Int64, posz);
        ny[1] = (ny[1] + Ncell)%Ncell+1; nz[1] = (nz[1] + Ncell)%Ncell+1;
    end
    
    nx, ny, nz = zeros(Int64, 2), zeros(Int64, 1), zeros(Int64, 1)
    wx, wy, wz = zeros(Float64, 2), zeros(Float64, 2), zeros(Float64, 2)
    @inbounds for i in 1:Nparticle
        mass = mass_p[i]
        FindxGwG_CIC!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        for j1 in 1:2
            CIC_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * mass
    end;end
    return CIC_array
end
# -----------------------------------------------------------------------------

#Triangle shape cloud, with equal particle mass.
#Interpolate the particle mass onto nearest three grids along each dim.
function TSC(posx_p, posy_p, posz_p, L, Ncell)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #end
    TSC_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function TSC_s(s)
        if abs(s) < 0.5
            return 3/4 - s^2
        elseif 0.5 <= abs(s) < 1.5
            return 0.5*(1.5-abs(s))^2
        else
            return 0
    end;end
    
    #Find the nearest three grid points around the particle along each dimension, with abs(s) < 1.5
    function FindxGwG_TSC!(posx::Float64, posy::Float64, posz::Float64, 
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1;
        s1 = n1 - posx;  s2 = n2 - posx;
        n3 = ifelse(s1 < -0.5, n2 + 1, n1 - 1); 
        s3 = n3 - posx
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1; nx[3] = (n3 + Ncell)%Ncell+1;
        wx[1] = TSC_s(s1);  wx[2] = TSC_s(s2);  wx[3] = TSC_s(s3);  
        
        n1 = floor(Int, posy); n2 = n1 + 1;
        s1 = n1 - posy;  s2 = n2 - posy;
        n3 = ifelse(s1 < -0.5, n2 + 1, n1 - 1); 
        s3 = n3 - posy
        ny[1] = (n1 + Ncell)%Ncell+1; ny[2] = (n2 + Ncell)%Ncell+1; ny[3] = (n3 + Ncell)%Ncell+1;
        wy[1] = TSC_s(s1);  wy[2] = TSC_s(s2);  wy[3] = TSC_s(s3);  
        
        n1 = floor(Int, posz); n2 = n1 + 1;
        s1 = n1 - posz;  s2 = n2 - posz;
        n3 = ifelse(s1 < -0.5, n2 + 1, n1 - 1); 
        s3 = n3 - posz
        nz[1] = (n1 + Ncell)%Ncell+1; nz[2] = (n2 + Ncell)%Ncell+1; nz[3] = (n3 + Ncell)%Ncell+1;
        wz[1] = TSC_s(s1);  wz[2] = TSC_s(s2);  wz[3] = TSC_s(s3);  
    end
    
    nx, ny, nz = zeros(Int64, 3), zeros(Int64, 3), zeros(Int64, 3)
    wx, wy, wz = zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3)
    @inbounds for i in 1:Nparticle
        FindxGwG_TSC!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx, wy, wz)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        for j1 in 1:3
            for j2 in 1:3
                for j3 in 1:3
                    TSC_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * wy[j2] * wz[j3]
    end;end;end;end
    return TSC_array
end

function TSC_w(posx_p, posy_p, posz_p, L, Ncell, mass_p)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #    masslist = reshape(mass_p, Nparticle);
    #end
    TSC_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function TSC_s(s)
        if abs(s) < 0.5
            return 3/4 - s^2
        elseif 0.5 <= abs(s) < 1.5
            return 0.5*(1.5-abs(s))^2
        else
            return 0
    end;end
    
    #Find the nearest three grid points around the particle along each dimension, with abs(s) < 1.5
    function FindxGwG_TSC!(posx::Float64, posy::Float64, posz::Float64, 
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1;
        s1 = n1 - posx;  s2 = n2 - posx;
        n3 = ifelse(s1 < -0.5, n2 + 1, n1 - 1); 
        s3 = n3 - posx
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1; nx[3] = (n3 + Ncell)%Ncell+1;
        wx[1] = TSC_s(s1);  wx[2] = TSC_s(s2);  wx[3] = TSC_s(s3);  
        
        n1 = floor(Int, posy); n2 = n1 + 1;
        s1 = n1 - posy;  s2 = n2 - posy;
        n3 = ifelse(s1 < -0.5, n2 + 1, n1 - 1); 
        s3 = n3 - posy
        ny[1] = (n1 + Ncell)%Ncell+1; ny[2] = (n2 + Ncell)%Ncell+1; ny[3] = (n3 + Ncell)%Ncell+1;
        wy[1] = TSC_s(s1);  wy[2] = TSC_s(s2);  wy[3] = TSC_s(s3);  
        
        n1 = floor(Int, posz); n2 = n1 + 1;
        s1 = n1 - posz;  s2 = n2 - posz;
        n3 = ifelse(s1 < -0.5, n2 + 1, n1 - 1); 
        s3 = n3 - posz
        nz[1] = (n1 + Ncell)%Ncell+1; nz[2] = (n2 + Ncell)%Ncell+1; nz[3] = (n3 + Ncell)%Ncell+1;
        wz[1] = TSC_s(s1);  wz[2] = TSC_s(s2);  wz[3] = TSC_s(s3);  
    end
    
    nx, ny, nz = zeros(Int64, 3), zeros(Int64, 3), zeros(Int64, 3)
    wx, wy, wz = zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3)
    @inbounds for i in 1:Nparticle
        mass = mass_p[i]
        FindxGwG_TSC!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx, wy, wz)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        for j1 in 1:3
            for j2 in 1:3
                for j3 in 1:3
                    TSC_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * wy[j2] * wz[j3] * mass
    end;end;end;end
    return TSC_array
end

function TSC1_w(posx_p, posy_p, posz_p, L, Ncell, mass_p)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #    masslist = reshape(mass_p, Nparticle);
    #end
    TSC_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function TSC_s(s)
        if abs(s) < 0.5
            return 3/4 - s^2
        elseif 0.5 <= abs(s) < 1.5
            return 0.5*(1.5-abs(s))^2
        else
            return 0
    end;end
    
    #Find the nearest three grid points around the particle along each dimension, with abs(s) < 1.5
    function FindxGwG_TSC!(posx::Float64, posy::Float64, posz::Float64, 
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1;
        s1 = n1 - posx;  s2 = n2 - posx;
        n3 = ifelse(s1 < -0.5, n2 + 1, n1 - 1); 
        s3 = n3 - posx
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1; nx[3] = (n3 + Ncell)%Ncell+1;
        wx[1] = TSC_s(s1);  wx[2] = TSC_s(s2);  wx[3] = TSC_s(s3);  
        
        ny[1] = round(Int64, posy);  nz[1] = round(Int64, posz);
        ny[1] = (ny[1] + Ncell)%Ncell+1; nz[1] = (nz[1] + Ncell)%Ncell+1;
    end
    
    nx, ny, nz = zeros(Int64, 3), zeros(Int64, 1), zeros(Int64, 1)
    wx, wy, wz = zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3)
    @inbounds for i in 1:Nparticle
        mass = mass_p[i]
        FindxGwG_TSC!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx, wy, wz)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        for j1 in 1:3
            for j2 in 1:3
                for j3 in 1:3
                    TSC_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * mass
    end;end;end;end
    return TSC_array
end
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



#Piecewise Cubic Spline, with equal particle mass.
#Interpolate the particle mass onto nearest four grids along each dim.
function PCS(posx_p, posy_p, posz_p, L, Ncell)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #end
    PCS_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function PCS_s(s)
        if s < 1             return (4 - 6s^2 + 3*s^3) /6
        elseif 1 <= s < 2      return (2 - s)^3 /6
        else                 return 0;      end
    end;
    
    #Find the nearest three grid points around the particle along a dimension, with abs(s) < 1.5
    function FindxGwG_PCS!(posx::Float64, posy::Float64, posz::Float64, 
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1; n3 = n1 + 2; n4 = n1 - 1;
        #Here s is positively defined
        s1 = posx - n1;  s2 = n2 - posx; s3 = n3 - posx;  s4 = posx- n4;
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1; nx[3] = (n3 + Ncell)%Ncell+1; nx[4] = (n4 + Ncell)%Ncell+1;
        wx[1] = PCS_s(s1);  wx[2] = PCS_s(s2);  wx[3] = PCS_s(s3);  wx[4] = PCS_s(s4)  
        
        n1 = floor(Int, posy); n2 = n1 + 1; n3 = n1 + 2; n4 = n1 - 1;
        #Here s is positively defined
        s1 = posy - n1;  s2 = n2 - posy; s3 = n3 - posy;  s4 = posy- n4;
        ny[1] = (n1 + Ncell)%Ncell+1; ny[2] = (n2 + Ncell)%Ncell+1; ny[3] = (n3 + Ncell)%Ncell+1; ny[4] = (n4 + Ncell)%Ncell+1;
        wy[1] = PCS_s(s1);  wy[2] = PCS_s(s2);  wy[3] = PCS_s(s3);  wy[4] = PCS_s(s4)  
      
        n1 = floor(Int, posz); n2 = n1 + 1; n3 = n1 + 2; n4 = n1 - 1;
        #Here s is positively defined
        s1 = posz - n1;  s2 = n2 - posz; s3 = n3 - posz;  s4 = posz- n4;
        nz[1] = (n1 + Ncell)%Ncell+1; nz[2] = (n2 + Ncell)%Ncell+1; nz[3] = (n3 + Ncell)%Ncell+1; nz[4] = (n4 + Ncell)%Ncell+1;
        wz[1] = PCS_s(s1);  wz[2] = PCS_s(s2);  wz[3] = PCS_s(s3);  wz[4] = PCS_s(s4)  
    end
    
    nx, ny, nz = zeros(Int64, 4), zeros(Int64, 4), zeros(Int64, 4)
    wx, wy, wz = zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4)
    @inbounds for i in 1:Nparticle
        FindxGwG_PCS!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx, wy, wz)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        @inbounds for j1 in 1:4
            @inbounds for j2 in 1:4
                @inbounds for j3 in 1:4
                    PCS_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * wy[j2] * wz[j3]
    end;end;end;end
    return PCS_array
end
# -----------------------------------------------------------------------------
#PCS with weight. Can be used to interpolate velocity or particles with different mass
function PCS_w(posx_p, posy_p, posz_p, L, Ncell, mass_p)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #    masslist = reshape(mass_p, Nparticle);
    #end
    PCS_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function PCS_s(s)
        if s < 1             return (4 - 6s^2 + 3*s^3) /6
        elseif 1 <= s < 2      return (2 - s)^3 /6
        else                 return 0;      end
    end;
    
    #Find the nearest three grid points around the particle along a dimension, with abs(s) < 2
    function FindxGwG_PCS!(posx::Float64, posy::Float64, posz::Float64,
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1}, wy::Array{Float64,1},  wz::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1; n3 = n1 + 2; n4 = n1 - 1;
        #Here s is positively defined
        s1 = posx - n1;  s2 = n2 - posx; s3 = n3 - posx;  s4 = posx- n4;
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1; nx[3] = (n3 + Ncell)%Ncell+1; nx[4] = (n4 + Ncell)%Ncell+1;
        wx[1] = PCS_s(s1);  wx[2] = PCS_s(s2);  wx[3] = PCS_s(s3);  wx[4] = PCS_s(s4)  
        
        n1 = floor(Int, posy); n2 = n1 + 1; n3 = n1 + 2; n4 = n1 - 1;
        #Here s is positively defined
        s1 = posy - n1;  s2 = n2 - posy; s3 = n3 - posy;  s4 = posy- n4;
        ny[1] = (n1 + Ncell)%Ncell+1; ny[2] = (n2 + Ncell)%Ncell+1; ny[3] = (n3 + Ncell)%Ncell+1; ny[4] = (n4 + Ncell)%Ncell+1;
        wy[1] = PCS_s(s1);  wy[2] = PCS_s(s2);  wy[3] = PCS_s(s3);  wy[4] = PCS_s(s4)  
      
        n1 = floor(Int, posz); n2 = n1 + 1; n3 = n1 + 2; n4 = n1 - 1;
        #Here s is positively defined
        s1 = posz - n1;  s2 = n2 - posz; s3 = n3 - posz;  s4 = posz- n4;
        nz[1] = (n1 + Ncell)%Ncell+1; nz[2] = (n2 + Ncell)%Ncell+1; nz[3] = (n3 + Ncell)%Ncell+1; nz[4] = (n4 + Ncell)%Ncell+1;
        wz[1] = PCS_s(s1);  wz[2] = PCS_s(s2);  wz[3] = PCS_s(s3);  wz[4] = PCS_s(s4)  
    end
    
    nx, ny, nz = zeros(Int64, 4), zeros(Int64, 4), zeros(Int64, 4)
    wx, wy, wz = zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4)
    @inbounds for i in 1:Nparticle
        mass = mass_p[i]
        FindxGwG_PCS!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx, wy, wz)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        @inbounds for j1 in 1:4
            @inbounds for j2 in 1:4
                @inbounds for j3 in 1:4
                    PCS_array[nx[j1], ny[j2], nz[j3]] += wx[j1] * wy[j2] * wz[j3] * mass
    end;end;end;end
    return PCS_array
end
# =============================================================================

# -----------------------------------------------------------------------------
#PCS with weight in 1-direction. Can be used to interpolate velocity or particles with different mass along 1 direction
function PCS1_w(posx_p, posy_p, posz_p, L, Ncell, mass_p)
    #dim_p = ndims(posx_p)
    Nparticle = length(posx_p)
    H = L / Ncell
    #if dim_p > 1
    #    posx_p = reshape(posx_p, Nparticle); 
    #    posy_p = reshape(posy_p, Nparticle); 
    #    posz_p = reshape(posz_p, Nparticle);
    #    masslist = reshape(mass_p, Nparticle);
    #end
    PCS_array = zeros(Float64, Ncell, Ncell, Ncell)
    
    #The interpolation function W(s)
    #s is the distance between the particle and the grid point
    function PCS_s(s)
        if s < 1             return (4 - 6s^2 + 3*s^3) /6
        elseif 1 <= s < 2      return (2 - s)^3 /6
        else                 return 0;      end
    end;
    
    #Find the nearest three grid points around the particle along a dimension, with abs(s) < 2
    function FindxGwG_PCS1!(posx::Float64, posy::Float64, posz::Float64,
                           nx::Array{Int64,1}, ny::Array{Int64,1}, nz::Array{Int64,1},
                           wx::Array{Float64,1})
        n1 = floor(Int, posx); n2 = n1 + 1; n3 = n1 + 2; n4 = n1 - 1;
        #Here s is positively defined
        s1 = posx - n1;  s2 = n2 - posx; s3 = n3 - posx;  s4 = posx- n4;
        nx[1] = (n1 + Ncell)%Ncell+1; nx[2] = (n2 + Ncell)%Ncell+1; nx[3] = (n3 + Ncell)%Ncell+1; nx[4] = (n4 + Ncell)%Ncell+1;
        wx[1] = PCS_s(s1);  wx[2] = PCS_s(s2);  wx[3] = PCS_s(s3);  wx[4] = PCS_s(s4)  
        
        ny[1] = round(Int64, posy);  nz[1] = round(Int64, posz);
        ny[1] = (ny[1] + Ncell)%Ncell+1; nz[1] = (nz[1] + Ncell)%Ncell+1;
    end
    
    nx, ny, nz = zeros(Int64, 4), zeros(Int64, 1), zeros(Int64, 1)
    wx = zeros(Float64, 4)
    @inbounds for i in 1:Nparticle
        mass = mass_p[i]
        FindxGwG_PCS1!(posx_p[i]/H, posy_p[i]/H, posz_p[i]/H, nx, ny, nz, wx)
        #println("$(posx_p[i]/H), $(posy_p[i]/H), $(posz_p[i]/H), $nx, $ny, $nz, $wx, $wy, $wz")
        @inbounds for j1 in 1:4
            PCS_array[nx[j1], ny[1], nz[1]] += wx[j1] * mass
    end;end;
    return PCS_array
end
# =============================================================================
