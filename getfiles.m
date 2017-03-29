%% Files aus gegebenen Ordnern sammeln
%%
%% Input: fo=Ordner als Integer
%%        str=Suchstring
%%
%% Output: files=Dateinamen als Strings in einer Structure 

function [files]=getfiles(str,fo)

root='/home/vonpapen/PhD/Cassini/';

if nargin<1
    str='10m*.wtf';
end
if nargin<2 || isempty(fo)==1
    tmp=dir(root);
    for i=3:length(tmp)
        fo(i-2)=strread(tmp(i).name,'%u');
    end
    clear tmp
end

k=1;

for i=1:length(fo)
    f=dir([root, num2str(fo(i),'%05i'), '/', str]);
    if ~isempty(f)
        for j=1:length(f)
            n{k}={[root, num2str(fo(i),'%05i'), '/', f(j).name]};
            k=k+1;
        end
    end
end
for i=1:length(n)
    files{i}=cell2mat(n{i});
%         fprintf('%s\n',files{i})
end