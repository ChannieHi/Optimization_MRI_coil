clc;
close;
clear all;

addpath('C://femm42/mfiles/')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       Overview   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 변수 설정: 코일 디자인 영역의 각 직사각형 영역에 대해 연속 변수를 사용하여 전류 소스의 존재 여부와 크기를 나타냄 (6*16)
% 변수와 목적 함수의 연결: 변수와 목적 함수를 연결하기 위해 가중치를 사용
% 가중치 설정 방법: 소스의 위치에 따른 직사각형 영역의 전자기 필드의 변동을 나타내는 지표(y축 방향 magnetic flux density의 표준편차와 x축 방향 flux density의 합) 구하여 음수 취하고 각 영역에서의 가중치 할당
% 목적 함수 정의: 가중치와 변수의 곱을 최소화 하는 것을 목표
% 최적화 알고리즘: Gradient descent algorithm 을 사용

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define optimization parameter
iter1 = 5;
iter2 = 15;
object_parameter = {}

%% Open FEMM42 and define problem
for k = 1:5
    for index1 = 0:iter1
        for index2 = 0:iter2

            std_y = {};
            flux_x = {};
    
            %Openning femm program
            openfemm;
            
            %Create magnetic problem
            newdocument(0);
            mi_probdef(0,'millimeters','axi');
            
            %Define node
            wx = 100;
            wy = 100;
            
            p1 = [0,0
                wx,0
                3*wx,0
                6*wx,0
                20*wx,0
                20*wx,20*wy
                0,20*wy
                0,wy
                wx,wy
                wx,0
                3*wx,0
                3*wx,8*wy
                6*wx,8*wy
                6*wx,0
                0,0
                0,wx];
            Drawpolygon(p1);
            mi_addsegment(p1(size(p1,1),1),p1(size(p1,1),2),p1(1,1),p1(1,2));
            
            %X-direction net
             p2 = zeros(11, 2); 
            
            for i = 1:11
                if mod(i,2) == 1
                    p2(i, 1) = 300 + (i+1)/2*50; 
                else 
                    p2(i,1) = p2(i-1);
                end
                
                if mod(i, 4) >= 2
                    p2(i, 2) = 800; 
                else
                    p2(i, 2) = 0;  
                end
            end
            Drawpolygon(p2);
            
            %y-direction net
            p3 = zeros(31, 2); 
            
            for i = 1:31
                if mod(i, 4) >= 2
                    p3(i, 1) = 600; 
                else
                    p3(i, 1) = 300;  
                end
                if mod(i,2) == 1
                    p3(i, 2) = (i+1)/2*50; 
                else 
                    p3(i,2) = p3(i-1,2);
                end
            end
            Drawpolygon(p3);
            
            % Add material air & coil label 
            mi_getmaterial('Air');
            SetAir(1000,1000,30);
            SetAir(50,50,30);
            mi_addmaterial('Coil',1,1,0,-k,0,0,0,1,0,0,0,0,0);

            %Set boundary conditions
            mi_addboundprop('boundary',0,0,0,0,0,0,0,0,0,0,0)
            mi_selectsegment(500,2000);
            mi_selectsegment(2000,1000);
            mi_setsegmentprop('boundary',0,0,0,0);
            

           
            
            % construction(index1,index2);
            construction();
            SetCoil(325 + 50*index1, 25 + 50*index2,'Coil');
    
            %Save file
            mi_saveas('MRI_matlab.fem')
            
            %Create mesh & Solve
            mi_createmesh();
            mi_analyze(0);
            mi_loadsolution();
    
            %load flux density
            flux_density = {};
            for i = 1:3
                for j = 1:3
                    flux_density{i,j} = mo_getb(25*i,25*j);
                end
            end
            temp_yflux = [];
            temp_xflux= [];
            for i = 1:3
                for j= 1:3
                temp_yflux(i,j) = flux_density{i,j}(1,2);
                temp_xflux(i,j) = flux_density{i,j}(1,1);
                end
            end
            std_y = std2(temp_yflux);
            flux_x = sum(temp_xflux,'all');
            temp_parameter(index1+1,index2+1) = - (abs(std_y) + flux_x);
            object_parameter{k,1} = temp_parameter;
            closefemm();
        end
    end
end
%% save weight
save("weight","object_parameter");

%% Gradient descent optimization for finding local Minimum
clc
clear all
close all

%Load parameter & Define object function 
load('weight.mat');
for i = 1:5
    object_parameter{i,1} = reshape(object_parameter{i,1},[96,1]);
end
temp_x = cell2mat(object_parameter);
X_norm = reshape(normalize(temp_x),[5,96]);
X = mean(X_norm);

%Set hyperparameters
learning_rate = 0.01;
num_iterations = 1000;

% Gradient descent optimization
weights = gradient_descent(X, learning_rate, num_iterations);

%Matrix processing
sorted_array = sort(weights);
num_intervals = 5;
interval_boundaries = linspace(min(sorted_array), max(sorted_array), num_intervals + 1);


for i = 1:numel(weights)
    for j = 1:num_intervals
        if weights(i) >= interval_boundaries(j) && weights(i) < interval_boundaries(j+1)
            weights(i) = j-1;
            break;
        end
    end
end

weights = reshape(weights,[6,16]);

%% Result
%Openning femm program
        openfemm;
        
        %Create magnetic problem
        newdocument(0);
        mi_probdef(0,'millimeters','axi');
        
        %Define node
        wx = 100;
        wy = 100;
        
        p1 = [0,0
            wx,0
            3*wx,0
            6*wx,0
            20*wx,0
            20*wx,20*wy
            0,20*wy
            0,wy
            wx,wy
            wx,0
            3*wx,0
            3*wx,8*wy
            6*wx,8*wy
            6*wx,0
            0,0
            0,wx];
        Drawpolygon(p1);
        mi_addsegment(p1(size(p1,1),1),p1(size(p1,1),2),p1(1,1),p1(1,2));
        
        %X-direction net
         p2 = zeros(11, 2); 
        
        for i = 1:11
            if mod(i,2) == 1
                p2(i, 1) = 300 + (i+1)/2*50; 
            else 
                p2(i,1) = p2(i-1);
            end
            
            if mod(i, 4) >= 2
                p2(i, 2) = 800; 
            else
                p2(i, 2) = 0;  
            end
        end
        Drawpolygon(p2);
        
        %y-direction net
        p3 = zeros(31, 2); 
        
        for i = 1:31
            if mod(i, 4) >= 2
                p3(i, 1) = 600; 
            else
                p3(i, 1) = 300;  
            end
            if mod(i,2) == 1
                p3(i, 2) = (i+1)/2*50; 
            else 
                p3(i,2) = p3(i-1,2);
            end
        end
        Drawpolygon(p3);
        
        % Add material air & coil label 
        mi_getmaterial('Air');
        SetAir(1000,1000,30);
        SetAir(50,50,30);
        mi_addmaterial('1',1,1,0,-0.1,0,0,0,1,0,0,0,0,0);
        mi_addmaterial('2',1,1,0,-0.2,0,0,0,1,0,0,0,0,0);
        mi_addmaterial('3',1,1,0,-0.3,0,0,0,1,0,0,0,0,0);
        mi_addmaterial('4',1,1,0,-0.4,0,0,0,1,0,0,0,0,0);

        %Set boundary conditions
        mi_addboundprop('boundary',0,0,0,0,0,0,0,0,0,0,0)
        mi_selectsegment(500,2000);
        mi_selectsegment(2000,1000);
        mi_setsegmentprop('boundary',0,0,0,0);
              
        % construction(index1,index2);
        construction();

        for i = 0:5
            for j = 0:15
                if weights(i+1,j+1) == 0
                   SetAir(325 + 50*i, 25+ 50*j,10);
                elseif weights(i+1,j+1) == 1
                   SetCoil(325 + 50*i,25 + 50*j,'1');
                   % SetAir(325 + 50*i, 25+ 50*j,10);
                elseif weights(i+1,j+1) == 2
                   SetCoil(325 + 50*i,25 + 50*j,'2');
                   % SetAir(325 + 50*i, 25+ 50*j,10);
                elseif weights(i+1,j+1) == 3
                   SetCoil(325 + 50*i,25 + 50*j,'3');
                else
                   SetCoil(325 + 50*i,25 + 50*j,'4');
                end
            end
        end
        %Save file
        mi_saveas('MRI_matlab.fem')
        
        %Create mesh & Solve
        mi_createmesh();
        mi_analyze(0);
        mi_loadsolution();

        %load flux density
        flux_density = {};
        for i = 1:3
            for j = 1:3
                flux_density{i,j} = mo_getb(25*i,25*j);
            end
        end
        temp_yflux = [];
        temp_xflux= [];
        for i = 1:3
            for j= 1:3
            temp_yflux(i,j) = flux_density{i,j}(1,2);
            temp_xflux(i,j) = flux_density{i,j}(1,1);
            end
        end
        std_y = std2(temp_yflux);
        flux_x = sum(temp_xflux,'all');
        result = (abs(std_y) + flux_x)

%% Functions
function Drawpolygon(p)
    for i = 1:size(p,1)
        mi_addnode(p(i,1),p(i,2));
    end
    for i = 1:size(p,1)-1
        mi_addsegment(p(i,1),p(i,2),p(i+1,1),p(i+1,2));
    end
    
end

function construction()
    for i = 0:5
        for j = 0:15
            SetAir(325 + 50*i, 25+ 50*j,10);
        end
    end
    
    % matrix = cell(11,31);
% for i = 1:11
%     for j = 1:31
%         if isempty(matrix{i,j})
%             matrix{i,j} = char(randi([111,112])); 
%             if strcmp(matrix(i,j),'o')
%             SetCoil(300 + 25*i,25*j)
%             end
%             if strcmp(matrix(i,j),'p')
%             SetAir(300+ 25*i,25*j,10)
%             end
%         end
%     end
% end
end

function SetAir(i,j,size)
    mi_addblocklabel(i,j);
    mi_selectlabel(i,j);
    mi_setblockprop('Air',0,size,'None',0,0,1);
    mi_clearselected();
end

function SetCoil(i,j,magnitude)
    mi_addblocklabel(i,j);
    mi_selectlabel(i,j);
    mi_setblockprop(magnitude,0,1,'None',0,0,1);
    mi_clearselected();
end

function weights = gradient_descent(X, learning_rate, num_iterations)
    num_features = size(X, 1);  
    weights = zeros(num_features, 1);

    for iter = 1:num_iterations
        predictions = X' .* weights;
        
        errors = predictions - X';
        
        gradient = 2 * X .* errors';
        
        weights = weights - learning_rate * gradient';
   end
end
