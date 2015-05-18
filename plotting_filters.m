function [filters,lpal] = plotting_filters(filters)

    if nargin < 2
		N = [];	
	end
	
	if isempty(N) && isfield(filters,'meta') && isfield(filters.meta,'size_filter')
		N = filters.meta.size_filter;
	else
		error('Unable to find max filter size!');
	end
	
	if length(N) == 1
		N = [N 1];
	end

	for p = 1:numel(filters.psi.filter)
		filter_coefft = realize_filter(filters.psi.filter{p},N);
		plot( abs(filter_coefft).^2/2 );hold on;
	end

	filter_coefft = realize_filter(filters.phi.filter,N);
	plot( abs(filter_coefft).^2/2 );hold on;
end 


