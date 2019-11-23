% read iterations

datafolder="~/lib/data/work-tiny-steps/work-tiny-steps-max_corr_steps_3-epsilon_0.0000001"
cd(datafolder);
myfiles=dir("*stderr*");
nfiles = size(myfiles,1);

niter_good = zeros(312*nfiles); % niter_good(i) = 1 if it corresponds to a good, physical solution %= zeros(1,nfiles);
niter = zeros(312*nfiles); %= zeros(311,nfiles);
for f=1:nfiles
  [~,niter_f] = system(['grep number\ of\ steps\ in\ path ' myfiles(f).name '|' 'cut -f 1 -d %'])
  eval (['niter(((f-1)*312+1):f*312) = [' niter_f ']']);
  % read found solution at index
end

% niter_good =[] %= zeros(1,nfiles);
% niter_bad  =[] %= zeros(311,nfiles);

% niter(nfiles,312)

% histogram niter_good
% histogram niter_bad
