function [ distance_matrix, pore_center ] = npc_get_opposite_su_distances(motl,bin_factor)
%% Compute distances of opposite subunits in the motl and pore center computed
% as an intersection of lines connecting opposite subunits (only if there
% are more than 1 line
%
%   Input:
%           motl: motive list specified as 20xN matrix
%           bin_factor: coordinates in motl will be scaled with this
%               number. If not given, value of 1 (i.e. no binning) is used
%
%   Output:
%           distance_matrix: 5xM matrix where M is number of pores in the
%               motl. First row contains pore number (motl row #6),
%               subsequent rows contain distances of opposite subunits (1-5, 
%               2-6, 3-7, and 4-8). If a pair of subunit does not exist,
%               the corresponding distance is set to 0.
%           pore_center: 8xM matrix where M is number of pores. First row
%               contains pore number, following 3 rows x,y,z, coordinates
%               of the center. Last 4 rows contain distances of respective lines
%               to the center. Center of a pore is computed as an intersection
%               of 3d lines connecting the opposite subunits. If there if
%               there is only 1 line (i.e. one pair of opposite subunits),
%               zeros are returned.
%

motl=emread(motl);
pores=unique(motl(6,:));

distance_matrix=zeros(5,numel(pores));
pore_center=zeros(8,numel(pores));

if(nargin==1)
    bin_factor=1;
end

for p=1:numel(pores)
    
    % take current pore motl
    npc_motl=motl(:,motl(6,:)==pores(p));
	
    %display(num2str(pores(p)));
    
	% sort npc motl according to su number
	npc_motl=sortrows(npc_motl',3)';
	
	% get subunit count
	su_count=size(npc_motl,2);
    
    distance_matrix(1,p)=pores(p);
    pore_center(1,p)=pores(p);
        
    if(su_count>=3)

        % take pore center coordinates
        centers=(npc_motl(8:10,:)+npc_motl(11:13,:)).*bin_factor;
    
        % create starting and ending points for lines defined by centers of
        % opposite subunits (1-5, 2-6, 3-7, 4-8)
        starting_points=[];
        ending_points=[];

        
        for j=1:min(4,su_count)

            su_idx1=j;
            su_idx2=find(npc_motl(3,:)==(npc_motl(3,j)+4));

            if(~isempty(su_idx2))

                starting_points = [ starting_points centers(:,su_idx1)];
                ending_points =[ ending_points centers(:,su_idx2)];

                % compute distances between opposite subunits
                distance_matrix(su_idx1+1,p)=sqrt(sum((centers(:,su_idx1)-centers(:,su_idx2)).^2));

            end


        end
        
        if(size(starting_points,2)>=2)
            
            % get the center of the lines 
            [center_fit, distances ] = geometry_3d_line_intersection(starting_points',ending_points');
            pore_center(2:4,p)=center_fit';
            
            dist_count=1;
            for d=2:5  
                if(distance_matrix(d,p)~=0)
                    pore_center(d+3,p)=distances(dist_count);
                    dist_count=dist_count+1;
                end
            end
            
        end
    
    end
    
end
		
end % function

function I = emread(motl_name)
%% Function was originaly part of TOM/AV3 package - was rewritten and slighlty 
% adapted to make above function stand-alone
% It reads motive list stored as EM format and specified by motl_name

% get file handle
fid = fopen(motl_name,'r','ieee-le');

% get file info
machine = fread(fid,[1],'uint8');
header = fread(fid,[2],'uint8');
data_type = fread(fid,[1],'uint8');

% close file
fclose(fid);

% check the machine and get endian based on that
if machine==6
    fid = fopen(motl_name,'r','ieee-le');
elseif machine==3
    fid = fopen(motl_name,'r','ieee-be');
elseif machine==5
    fid = fopen(motl_name,'r','ieee-be');
else
    error('Error: Wrong File Format');
end

header = fread(fid,[128],'uint32');
xdim = header(2);
ydim = header(3);
zdim = header(4);

if data_type==5
    I = reshape(fread(fid,[xdim * ydim * zdim],'float'),xdim, ydim, zdim);
elseif data_type==1
    I = uint8(reshape(fread(fid,[xdim * ydim * zdim],'uint8'),xdim, ydim, zdim));
elseif data_type==2
    I = uint8(reshape(fread(fid,[xdim * ydim * zdim],'uint16'),xdim, ydim, zdim));
else
    error('Error: Wrong Data Type');
end

fclose(fid);

end
