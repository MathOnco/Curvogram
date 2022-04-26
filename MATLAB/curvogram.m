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
%           'Color' - bin color
%           'FaceAlpha' - alpha value (opacity) of bins
%           'SkewBins' - true or false (false = straight, parallel lines)
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
BinWidth = 1;
ScaleHeight=1;
Color = [0.2188,0.4531,0.6914]; % default is blue
FaceAlpha = 1.0;

% set a buffer on default xlimits:
xlimits = [min(X) - (max(X)-min(X))*0.5, max(X) + (max(X)-min(X))*0.5];

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

SkewBins=true;


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
    elseif strcmp(varargin{index},'SkewBins')
        SkewBins=varargin{index+1};
    elseif strcmp(varargin{index},'Color')
        Color=varargin{index+1};
    elseif strcmp(varargin{index},'FaceAlpha')
        FaceAlpha=varargin{index+1};
    end
end


% add gridlines at same locations as the histogram:
addGridlines(X,Y,f,fprime,xlimits,ylimits,ScaleHeight);

% plot function for full range:
plot(x,y,'-k','LineWidth',5); hold on;

custom_histogram(f,fprime,X,Y,xlimits,ylimits,BinWidth,ScaleHeight,SkewBins,Color,FaceAlpha);


xlim(xlimits);
ylim(ylimits);

end


%% add gridlines
function [] = addGridlines(X,Y,f,fprime,xlimits,ylimits,ScaleHeight)

    SCALE = diff(xlimits)/diff(ylimits);

    d = (max(Y)*ScaleHeight)*1.1;    

    % place gridlines on "X" vals, plus more on either side:
    xStep = mean(min(diff(X),max(diff(X)))); % average step size of "X"
    xVecCoarse = [X(1)-xStep,X,X(end)+xStep];

    % assume 10 grid lines:
    xStepSmall = (xVecCoarse(end)-xVecCoarse(1))/(200);
    xVecFine = xVecCoarse(1): xStepSmall : xVecCoarse(end);
    yVec = 0 : (d/10) : d;


    gray = [227,227,227]/350;

    %% add vertical lines:
    for x = xVecCoarse
        point=[x,d]; % coordinates in transformed system
        P=coord_transform(point,SCALE,f,fprime);
        plot([x,P(1)],[f(x),P(2)],'-','LineWidth',1,'Color',gray); hold on;
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
        plot(line_x,line_y,'-','LineWidth',1,'Color',gray); hold on;
    
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



%% add polygons at each y location
function [] = custom_histogram(f,fprime,X,Y,xlimits,ylimits,BinWidth,ScaleHeight,SkewBins,Color,FaceAlpha)

    % colors:
    red = [0.84,0.18,0.14];
    blue = [0.22,0.45,0.69];
    
    % helper values
    SCALE = diff(xlimits)/diff(ylimits);
    margin = 1-BinWidth;% default is 0.1;
    delta_x = X(2) - X(1);
    w = delta_x*(1-margin);


    i=1;
    for xVal = X

        h=Y(i)*ScaleHeight;
        theta1 = atan2(fprime(xVal), 1/SCALE);        

        %% plot line only:
        bottom = coord_transform([xVal,0],SCALE,f,fprime); % bottom
        top = coord_transform([xVal,h],SCALE,f,fprime); % top

        % find corners if SkewBins is true:
        
        % bottom left:
        BL = coord_transform([xVal-w/2,0],SCALE,f,fprime);
        x2 = BL(1);
        y2 = BL(2);

        % bottom right:
        BR = coord_transform([xVal+w/2,0],SCALE,f,fprime);
        x3 = BR(1);
        y3 = BR(2);

        % top right:
        TR = coord_transform([xVal+w/2,h],SCALE,f,fprime);
        x4 = TR(1);
        y4 = TR(2);

        % top left:
        TL = coord_transform([xVal-w/2,h],SCALE,f,fprime);
        x1 = TL(1);
        y1 = TL(2);

        if (~SkewBins)
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
        end
        
        if (h>0)    
            pgon = polyshape([x1,x2,x3,x4],[y1,y2,y3,y4]);hold on;
            c = (red-blue).*(h/(max(Y)*ScaleHeight)) + blue;
            if ~isnan(Color)
                c=Color;
            end
            plot(pgon,'FaceAlpha',FaceAlpha,'FaceColor',c,'LineWidth',1.5);
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







