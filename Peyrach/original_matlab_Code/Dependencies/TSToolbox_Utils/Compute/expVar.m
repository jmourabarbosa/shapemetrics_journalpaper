function EV = expVar(r_w,r_pre,r_post)


% =========================================================================
%                            expVar
% =========================================================================
% 
% USAGE: EV = expVar(r_w,r_pre,r_post)
%
% DESCRIPTION:  Computes explained variance
%
% =========================================================================
% INPUTS: 
%    __________________________________________________________________
%       Properties          Description                     Default
%    __________________________________________________________________
%
%       r_w                 correlation coeff. from wake
%       r_pre               same for sleep pre
%       r_post              same for sleep post
%
% Adrien Peyrache

% Make sure we're dealing with vertical vectors here
r_w = r_w(:);
r_pre = r_pre(:);
r_post = r_post(:);

r_w_post    = nancorrcoef(r_w,r_post);
r_w_post    = r_w_post(1,2);
r_w_pre     = nancorrcoef(r_w,r_pre);
r_w_pre     = r_w_pre(1,2);
r_pre_post  = nancorrcoef(r_post,r_pre);
r_pre_post  = r_pre_post(1,2);

EV = (r_w_post - r_w_pre * r_pre_post);
EV = EV / sqrt((1 - r_w_pre^2)*(1 - r_pre_post^2));
EV = EV ^ 2;
