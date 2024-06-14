using Base
using HomotopyContinuation
using LinearAlgebra

#Rotation parametrization
function cay2R(x,y,z)
    M = [1+x*x-(y*y+z*z)  2*(x*y-z)  2*(x*z+y);
         2*(x*y+z)  1+y^2-(x*x+z*z)  2*(y*z-x);
         2*(x*z-y)  2*(y*z+x)  1+z*z-(x*x+y*y)];
    return M;
end 
@var r2[1:3] r3[1:3]
R2 = cay2R(r2[1], r2[2], r2[3]);
R3 = cay2R(r3[1], r3[2], r3[3]);

#Set translations
@var t2[1:3] t3[1:3]

#Set parameters
@var x[1:3,1:3,1:3];
@var d[1:2,1:3,1:3];

#Set unknowns
@var a[1:3,1:3];
@var e[1:2,1:3];
@var u[1:2,1:3];

#charts
@var c[1:6];

#Point Equations
p = 1
pointEquations2 = a[p,2] * x[p,2,:] - (R2 * (a[p,1] * x[p,1,:]) + t2)
for p = 2:3
    Eq = a[p,2] * x[p,2,:] - (R2 * (a[p,1] * x[p,1,:]) + t2)
    for n = 1:3
        push!(pointEquations2, Eq[n]);
    end
end
p = 1
pointEquations3 = a[p,3] * x[p,3,:] - (R3 * (a[p,1] * x[p,1,:]) + t3)
for p = 2:3
    Eq = a[p,3] * x[p,3,:] - (R3 * (a[p,1] * x[p,1,:]) + t3)
    for n = 1:3
        push!(pointEquations3, Eq[n]);
    end
end
pointEquations = [pointEquations2; pointEquations3];

#Tangent Equations
p = 1
tangentEquations2 = (e[p,2] * x[p,2,:] + u[p,2] * d[p,2,:]) - (R2 * (e[p,1] * x[p,1,:] + u[p,1] * d[p,1,:]));
for p = 2:2
    Eq = (e[p,2] * x[p,2,:] + u[p,2] * d[p,2,:]) - (R2 * (e[p,1] * x[p,1,:] + u[p,1] * d[p,1,:]))
    for n = 1:3
        push!(tangentEquations2, Eq[n]);
    end
end
p = 1
tangentEquations3 = (e[p,3] * x[p,3,:] + u[p,3] * d[p,3,:]) - (R3 * (e[p,1] * x[p,1,:] + u[p,1] * d[p,1,:]));
for p = 2:2
    Eq = (e[p,3] * x[p,3,:] + u[p,3] * d[p,3,:]) - (R3 * (e[p,1] * x[p,1,:] + u[p,1] * d[p,1,:]))
    for n = 1:3
        push!(tangentEquations3, Eq[n]);
    end
end
tangentEquations = [tangentEquations2; tangentEquations3];

#Charts Equations
charts = [c[1] * a[1,1] - c[2];
          c[3] * e[1,1] - c[4]; 
          c[5] * e[2,1] - c[6]];

#All Equations, 33 equations in total
Eqs = [pointEquations; tangentEquations; charts];
variables_list = collect(Iterators.flatten([transpose(a),transpose(e),transpose(u),t2,t3,r2,r3]));
parameters_list = collect(Iterators.flatten([permutedims(x,[3,2,1]),permutedims(d,[3,2,1]),c]))
F = System(Eqs;variables=variables_list, parameters =parameters_list);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

#Try to solve with random Parameters
p_rand = rand(51) + rand(51) * im

@time for i = 1:100
solve(F, start_solutions; start_parameters=start_params, target_parameters=p_rand)
end


# #Try to solve with random Parameters
# p_rand = rand(51) + rand(51) * im
# F_prand = subs(Eqs,parameters_list => p_rand);
# sol_xrand = solve(F_prand)





# #Try to Solve with generated seeds from Macaulay2 with the command "(p0, x0) = createSeedPair G"
# p0 = [-.106919-.00213814im,-.0943637+.0316458im,-.253802-.115537im,1.15814-.0209775im,1.58467-.0902769im,.705916+.186359im,.330181+.0443676im,-.359701-.149705im,.151317-.247041im,-.360448+.220939im,.101598+.0386087im,-.161494-.168534im,.284865-.344842im,-.0827812+.153128im,.0563866-.482724im,.465327-.503302im,-1.59768+.476787im,.423075-.46709im,-.0855705+.0160487im,-.0495171-.264141im,-.221294+.105942im,.870365+.228331im,.288691-.497529im,.129296-.716875im,-.157768+.82508im,-.603055-.168421im,.579734-1.03396im,-.0311147+.0120491im,.0923877-.0303556im,.221007+.0958708im,-.569672-.195251im,-.983581-.603664im,-.136429-.38835im,.898365+.0840939im,.225448+.884679im,.137016+.38008im,.569659+.280488im,.181062-.165189im,.223627+.291214im,-.112841+.795554im,-.249752+.436785im,.188865+.50069im,.262957+.597315im,.844966+.494435im,.821847+.560416im,.808464+.699031im,.46707+.661389im,.751717+.544675im,-.0198735+.334826im,.986081+.175158im,.640182+.945381im]
# x0 = [.735333+.182282im,.43882+.405319im,.52235+.190799im,.962496+.614761im,.0890666+.715857im,.682343+.77288im,.855443+.06981im,.140385+.959576im,.961344+.167995im,.194293+.304635im,.313954+.467125im,.409802+.630922im,.794451+.817608im,.416113+.989015im,.554021+.316676im,.549325+.413572im,.931278+.475234im,.0306418+.0649098im,.570449+.248149im,.0140022+.973551im,.264005+.263418im,.52948+.932798im,.699406+.269814im,.526089+.591768im,.526912+.813107im,.0226996+.0234142im,.34661+.0717123im,.971845+.557098im,.273048+.015302im,.666807+.510207im,.97162+.345702im,.986023+.135618im,.935078+.44708im]
# F_p0 = subs(Eqs,parameters_list => p0);
# sol_x0 = solve(F_p0);

# counter = 0;
# flag = 0;
# for i = 1:376
#     global counter;
#     global flag;
#     flag = 0;
#     for j = 1:33
#         global counter;
#         global flag;
#         N = norm(sols_x0[i][j]);
#         if N < 1e-10
#             flag = 10;
#         end
#     end
#     if flag == 10
#         counter+=1;
#     end
#     flag = 0;
# end




