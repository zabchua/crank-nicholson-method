// CRANK-NICHOLSON ALGORITHM
// By Mary Elizabeth E. Chua
// Plug-and-play so change the length, dx, Tl, Tr, dt, and time. If the whole rod doesn't start at 0 then configure the U() array.
// The constants C, p, and k' (k1) is changeable as well

// How long is the rod?
length = 10

// What is delta x?
dx = 0.5

// Then the number of nodes is
m = (length/dx) - 1

// For the matrix we have n equations including the bounds
n = m + 2

// The two bounds, btw are
Tl = 100
Tr = 50

// The time interval shall be
dt = 0.01

// Until
time = 10
count = time/dt
count1 = count

// Then the starting temp for the interior nodes are
U = zeros(count+1, m+2)
U(:, 1) = Tl
U(:, m+2) = Tr
//U() = something if the rod doesn't start at 0

// Furthermore, here are the constants for the material
C = 0.835
p = 2.7
k1 = 0.49

// Therefore, k is
divisor = C*p
k = k1/divisor

// And lambda is
lambda = (k*dt)/(dx^2)

T = zeros(n, n+1);
double(T)

T(1, 1) = 1
T(1, n+1) = Tl
T(n, n) = 1
T(n, n+1) = Tr

while(count~=0)
    for (j=2:n-1)
        T(j, j-1) = -lambda
        T(j, j) = 2*(1+lambda)
        T(j, j+1) = -lambda
        
        if(j==2)
            T(j, n+1) = lambda*U(count+1,1) + 2*(1-lambda)*U(count+1, 2) + lambda*U(count+1, 3) 
        elseif(j==n-1)
            T(j, n+1) = lambda*U(count+1, n) + 2*(1-lambda)*U(count+1, n-1) + lambda*U(count+1, n-2) 
        else
            T(j, n+1) = lambda*U(count+1, j-1) + 2*(1-lambda)*U(count+1, j) + lambda*U(count+1, j-1)
        end
    end
    
//    disp(T)
    solve = rref(T)
//    disp(solve)
    
    for(i=2:m+1)
        U(count, i) = solve(i, n+1)
    end
    count = count-1
end

disp(U)

//ANIMATE THE GRAPH OVER TIME
//Only for problem 7
//help = U
//fprintfMat(TMPDIR + "/Mat", help, "%f");
//clear
//a1 = fscanfMat(TMPDIR + '/Mat');
//hat = 0:1000
//dog = 0:20
//clf()
//gcf().color_map = hotcolormap(256);
//
//for (i=1000:-1:0)
//    clf()
//    interval = (i/100)
//    Sgrayplot(hat, dog, a1, strf="011", rect=[0,i/100,10,interval+0.01])
//    sleep(10)
//end
