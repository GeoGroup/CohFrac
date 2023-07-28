function out = DFN(varargin)
% DFN
% generates fracture network models, 2D|3D
%
% Usage...:
% out = DFN(varargin);
%
% Input...: varargin  any
% Output..: out       struct,{Line,Orig|Poly,Orig}
%
% Examples:
%{
fnm = DFN; % 2d fnm
fnm = DFN('dim',2,'n',300,'dir',45,'ddir',0,'minl',0.1,'mu',0.2,'maxl',0.3); % 2d fnm
fnm = [Field(DFN('dim',2,'n',30,'dir',45,'ddir',-100,'minl',0.1,'mu',0.2,... % combined fnm
    'maxl',0.3),'Line');Field(DFN('dim',2,'n',30,'dir',135,'ddir',-100,'minl',...
    0.1,'mu',0.2,'maxl',0.3),'Line')];
fnm = Field(DFN('dim',3),'Poly'); % 3d fnm
fnm = Field(DFN('dim',3,'dip',45,'ddip',0,'dir',180,'ddir',0),'Poly'); % 3d fnm
%}
%
% Alghalandis Discrete Fracture Network Engineering (ADFNE),*R1.5*
% Copyright (c) 2018 Alghalandis Computing @
% Author: Dr. Younes Fadakar Alghalandis
% (w) http://alghalandis.net        (e) alghalandis.net@gmail.com
% All rights reserved.
%
% License.: ADFNE1.5_License.txt and at http://alghalandis.net/products/adfne/adfne15
%
% Citations:
% Fadakar-A Y, 2017, "ADFNE: Open source software for discrete fracture network
% engineering, two and three dimensional applications", Journal of Computers &
% Geosciences, 102:1-11.
%
% Fadakar-A Y, 2018, "DFNE Practices with ADFNE", Alghalandis Computing, Toronto, 
% Ontario, Canada, http://alghalandis.net, pp61.
%
% see more at: http://alghalandis.net/products/adfne
% Updated.: 2018-01-11

global Report
opt = Option(varargin,'n',100,'minl',0.05,'mu',0.3,'maxl',0.6,...               % default arguments
    'bbx',[0,0,1,1],'dim',2,'asep',0,'dsep',0,'mit',100,'scale',1,...           %n是线段的数量
    'shape','c','q',24,'dip',45,'ddip',-1e-7,'dir',0,'ddir',-1e-7);
Ticot(sprintf('Generating %dD FNM [n = %d]',opt.dim,opt.n));                    % initializes timing
lhs = 0.5*Rand(opt.n,'fun','exp','mu',opt.mu,'ab',[opt.minl,opt.maxl]); % half lengths~ Exp(mu) dist.
new_ols=[];
if opt.dim == 2                                                                 % 2D case                                                                     % simple simulation
    n = opt.n;
    pts = rand(n,2);                                                        % locations~ U(0,1)
    for i=1:n-1
        while norm(pts(i,:)-pts(i+1,:)) < 0.2
            pts(i,:) = rand(1,2);
        end
    end
    ags = Rand(n,'fun','f','mu',opt.dir,'k',-opt.ddir);                     % orientation, Fisher dist. 
    [dx,dy] = pol2cart(Convert('rad',ags),lhs);
    ols = [pts(:,1)-dx,pts(:,2)-dy,pts(:,1)+dx,pts(:,2)+dy];                % original lines
    for i = 1:size(ols,1)
        juli=[];
        for j = 1:size(ols,1)
            dis = line_distance(ols(i,:),ols(j,:));
            if dis ~=0
                juli=[juli,dis];
            end
        end
        min_dis = min(juli);
        if min_dis > 0.008
            min_dis
            new_ols = [new_ols;ols(i,:)];
        end
    end
    out.Orig = new_ols;
    out.Line = Clip(new_ols,opt.bbx);                                               % clipped by region of study
    if Report                                                                   
        Display('Original Lines#',size(out.Orig,1),...                          % displays information
            'Clipped Lines#',size(out.Line,1));
    end
end
Ticot;                                                                          % ends timing
