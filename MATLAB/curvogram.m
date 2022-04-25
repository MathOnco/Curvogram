function [] = curvogram(X,Y,f,fprime,varargin)
%curvogram(X,Y,f,fprime)
%   Takes as input:
%       X - vector (single row) of x-values (histogram bins)
%       Y - vector (single row) of y-values (histogram counts)
%       f - function, f(x):
%           e.g. f = @(x) x.^2
%       fprime - derivative of f(x):
%           e.g. fprime = @(x) 2*x
%
%       optional arguments: 
%           'XLimits' - x axis limits ( e.g. [0, 1] )
%           'YLimits' - y axis limits ( e.g. [0, 1] )
%           'BinWidth' - usually between 0 and 1 (default:0.9)
%           'ScaleHeight' - 1 is default
%       
%       example:
%           X = ......
%           curvogram(X,Y,f,fprime) .........
%           curvogram([0,1,2,3,4],[2,3,4,5,6],1,2,'XLimits',[0 10])

% get current figure, if exists:
h = gcf;
figure_number=h.Number;
figure(figure_number); hold on;




%% assert length(X) == length(Y):
errorMsg1 = 'X and Y must be the same length.'; 
XYLength = @(x,y) assert(length(x)==length(y),errorMsg1);
XYLength(X,Y);


%% get optional input:
p = inputParser; 

%% set default values:
BinWidth = 0.9;
ScaleHeight=1;

% set a buffer on default xlimits:
xlimits = [min(X) - (max(X)-min(X))*0.25, max(X) + (max(X)-min(X))*0.25];

% function values:
x= xlimits(1):((xlimits(2)-xlimits(1))/100):xlimits(end);
y=f(x);

ylimits = [min(y), 2.1*(max(y)-min(y)) + min(y)]; % default is to make it twice the span (leave room for the curvogram)

if (diff(ylimits)==0)
    ylimits = [ylimits(1),ylimits(1)+1];
end

% validation of user input labels:
addParameter(p,'XLimits',xlimits);
addParameter(p,'YLimits',ylimits);



% read in optional parameters
[nParams] = length(varargin);
for param = 1:1:(nParams/2)
    index = (param-1)*2 + 1;
    if strcmp(varargin{index}, 'XLimits')
        xlimits=varargin{index+1};
    elseif strcmp(varargin{index}, 'YLimits')
        ylimits=varargin{index+1};
    elseif strcmp(varargin{index}, 'ScaleHeight')
        ScaleHeight=varargin{index+1};
    elseif strcmp(varargin{index}, 'BinWidth')
        BinWidth=varargin{index+1};
    end
end





% add gridlines at same locations as the histogram:
addGridlines(X,f,fprime,xlimits,ylimits);



% plot function for full range:
plot(x,y,'-k','LineWidth',5); hold on;

custom_histogram(f,fprime,X,Y,xlimits,ylimits,BinWidth,ScaleHeight);


xlim(xlimits);
ylim(ylimits);

end


%% add gridlines
function [] = addGridlines(X,f,fprime,xlimits,ylimits)

    SCALE = diff(xlimits)/diff(ylimits);

    % d:= max distance of any line w/ axis scale limits:
    d = sqrt((xlimits(2) - xlimits(1)).^2 + (ylimits(2) - ylimits(1)).^2);


    % assume 20 grid lines:
    GL = 20;
    xStepSmall = (xlimits(2)-xlimits(1))/(GL*5);

    xVecFine = xlimits(1): xStepSmall : xlimits(2);
    yVec = 0 : ((ylimits(2))/GL) : d;

    % place gridlines on "X" vals, plus more on either side:
    xStep = mean(min(diff(X),max(diff(X)))); % average step size of "X"
    x_left = X(1):-xStep:xlimits(1);
    x_left = fliplr(x_left);
    x_right = X(end):xStep:xlimits(2);

    xVecCoarse = [x_left(1:end-1),X,x_right(2:end)];


    %% add vertical lines:
    for x = xVecCoarse
        point=[x,d]; % coordinates in transformed system
        P=coord_transform(point,SCALE,f,fprime);
        plot([x,P(1)],[f(x),P(2)],'-k','LineWidth',1); hold on;
    end
    
    % spacing of horizontal lines:
    for y = yVec
    
        line_x=[];
        line_y=[];
    

        % distance along curve:
        for x = xVecFine
    
            P2=coord_transform([x,y],SCALE,f,fprime);
    
            if (~isnan(P2(1)))
                line_x=[line_x,P2(1)];
                line_y=[line_y,P2(2)];
            end 
        end
        plot(line_x,line_y,'-k','LineWidth',1); hold on;
    
    end


    

    xlim(xlimits);
    ylim(ylimits);
    

    % beautify the plot & make sure that it is square in dimensions
    % (required for bars to visually be perpendicular)
    fs = 20;
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fs,'FontWeight','Bold', 'LineWidth', 2);
    
    p = get(gcf,'Position');
    plot_size = max(p(3:4));
    set(gcf,'Position',[p(1),p(2),plot_size,plot_size]);
    box on;
    set(gcf,'color','w');

end




function [] = custom_histogram(f,fprime,x_values,y_values,xlimits,ylimits,BinWidth,ScaleHeight)

    % colors:
    red = [0.84,0.18,0.14];
    blue = [0.22,0.45,0.69];
    
    % helper values
    SCALE = diff(xlimits)/diff(ylimits);
    margin = 1-BinWidth;% default is 0.1;
    delta_x = x_values(2) - x_values(1);
    w = delta_x*(1-margin);


    i=1;
    for xVal = x_values

        h=y_values(i)*ScaleHeight;
        theta1 = atan2(fprime(xVal), 1/SCALE);        

        %% plot line only:
        bottom = coord_transform([xVal,0],SCALE,f,fprime); % bottom
        top = coord_transform([xVal,h],SCALE,f,fprime); % top


        % bottom left:
        x2 = bottom(1) - (w/2)*cos(theta1);
        y2 = bottom(2) - (1/SCALE)*(w/2)*sin(theta1);

        % bottom right:
        x3 = bottom(1) + (w/2)*cos(theta1);
        y3 = bottom(2) + (1/SCALE)*(w/2)*sin(theta1);

        % top right:
        x4 = top(1) + (w/2)*cos(theta1);
        y4 = top(2) + (1/SCALE)*(w/2)*sin(theta1);

        % top left:
        x1 = top(1) - (w/2)*cos(theta1);
        y1 = top(2) - (1/SCALE)*(w/2)*sin(theta1);
        
        if (h>0)    
            pgon = polyshape([x1,x2,x3,x4],[y1,y2,y3,y4]);hold on;
            c = (red-blue).*(h/(max(y_values)*ScaleHeight)) + blue;
            plot(pgon,'FaceAlpha',0.8,'FaceColor',c);
        end

        % plot function for full range:
        x= xlimits(1):((xlimits(2)-xlimits(1))/100):xlimits(end);
        plot(x,f(x),'-k','LineWidth',5); hold on;

        i=i+1;
    end
end


% p is the x coordinate in the new system (along the line)
% q is the y coordinate in the new system
function [new_point] = coord_transform(point,SCALE,f,fprime)
    
    % coordinate lengths:
    p = point(1); % this is actually same x
    q = point(2); % distance along perpendicular bisector
        
    % angle created by the slope of the line:
    theta1 = atan2(fprime(p), 1/SCALE);

    % transformed coordinates:
    new_point = [p - SCALE*q*sin(theta1), f(p) + q*cos(theta1)];
    
end







