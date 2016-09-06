%This programm is SelectThickness
%It chooses the columns which contain thicknesses (notice: the thickness is the third row of each column, 
%except for the first four columns) within the boundaries of the chosen thicknesses and
%saves them in a new file named DRS_ML_cut

function index=SelectThickness(DRS_ML, dstart, dstop);  %one parameter of the function is DRS_ML so it hasn't to be loaded again
thickness=DRS_ML(3,:);
thickness(1:2)=0;
[~,mini]=min(abs(thickness-dstart));
[~,maxi]=min(abs(thickness-dstop));
index.min=min(mini, maxi);
index.max=max(mini, maxi);
end
