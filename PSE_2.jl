using LinearAlgebra
using Plots

### define constants
K = 0.025
D = 0.063
a = 4.2
Kb = 8.617*10^-5

function V(x)
   D*(1-exp(-x))^2
end
function W(x,y)
 (K/(2*a^2))*(x-y)^2
end
function KER(x,y,tempe)
exp(-(  V(x)*0.5 + V(y)*0.5+ W(x,y) )/(Kb*tempe) )
end

#### create grid
ndiv  = 8000
lmin = -40
lmax = 5000
temp_discret=60
grid = range(lmin, stop = lmax, length=ndiv+1)
delta = (lmax-lmin)/ndiv

### create transfer matrix


function create_transf_fava(tempe)
  transfer=zeros(length(grid),length(grid))
    for i = 1 : length(grid)
	   for j = 1 : length(grid)
	     transfer[i,j]=delta^2*KER(grid[i],grid[j],tempe)
       end
	end
return transfer	
end	


tempe = range(1,stop=600,length=temp_discret)
bet = zeros(temp_discret)
for j = 1 : length(tempe)
   bet[j]=1/(Kb*tempe[j])
end

freeenergie = zeros(temp_discret)
w = zeros(temp_discret)
for k =1:length(tempe)
    transfer = zeros(length(grid),length(grid))
	for i = 1 : length(grid)
	   for j = 1 : length(grid)
	     transfer[i,j]=delta^2*KER(grid[i],grid[j],tempe[k])
       end
	end   
	w[k] = maximum(eigvals(transfer))
    freeenergie[k] =-(Kb*tempe[k])*log( w[k])
end


#### Algorithm #2

function create_transf_one(tempe)
  transfer_one=zeros(length(grid),length(grid))
    for i = 1 : length(grid)
	   for j = 1 : length(grid)
	     transfer_one[i,j]=KER(grid[i],grid[j],1)
       end
	end
return transfer_one	
end	

#freeenergie_one = zeros(temp_discret)
#eigen_one = maximum(eigvals(create_transf_one(1)))
#for k = 1:length(tempe)
#    freeenergie_one[k] =-(Kb*tempe[k])*(2*log(delta) + (1/tempe[k])*log(eigen_one) )
#end

##### Plotting
 	      
plot(tempe,freeenergie)	

#logw=zeros(temp_discret)
#for i = 1 : length(tempe)
#    logw[i] = log(w[i])
#end

#plot(tempe, logw)
	

