% read iterations
clear

% datafolder="~/lib/data/work-tiny-steps/work-tiny-steps-max_corr_steps_3-epsilon_0.0000001";
datafolder="~/lib/data/work-ccv-steps/work-100_triplets-no_gammify_start-steps-ccv-max_corr_steps_4-epsilon_0.0000001";
cd(datafolder);
myfiles=dir("*stderr*");
nfiles = size(myfiles,1);

niter_good = zeros(312,nfiles); % niter_good(i) = 1 if it corresponds to a good, physical solution %= zeros(1,nfiles);
niter = zeros(312,nfiles); %= zeros(311,nfiles);
for f=1:nfiles
  [~,niter_f] = system(['grep number\ of\ steps\ in\ path ' myfiles(f).name '|' 'cut -f 1 -d %']);
  eval (['niter(:,f) = [' niter_f '];']);
  
  [found,good_id_f] = system(['grep solution\ at\ index: ' myfiles(f).name '|' 'cut -f 2 -d :']);
  [foundtmp,good_id_tmp] = system(['grep solution\ at\ index: ' myfiles(f).name ])
  if (found == 0)
    good_id_f
    eval (['good_id = [' good_id_f ']']);
    niter_good(good_id+1,f) = 1;
  end
  % read found solution at index
end

ngl = niter_good(:);
nil = niter(:);

ng = nil(find(ngl));
nb = nil(~ngl);

% histogram niter_good
% histogram niter_bad
