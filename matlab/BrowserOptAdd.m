function [ Opt ] = BrowserOptAdd( Opt0, foldername, indexname )
%BROWSEROPTADD Summary of this function goes here
%   Detailed explanation goes here
Opt = Opt0;
F = cell(1,1);
I = cell(1,1);
F{1} = foldername;
I{1} = indexname;

Opt.n = Opt.n+1;
Opt.index = [Opt.index;I];
Opt.folder = [Opt.folder;F];


end

