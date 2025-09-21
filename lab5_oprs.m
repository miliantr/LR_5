close all
maximum = 5;

%Project shadow

%r = read_data("res1.txt",4)

% l = size(r)
% for i = [1:l(1)]
%     plot3(r(i,3),r(i,1),r(i,2),'-r.');
%     h = r(i,5)/3600.0;
%     angle = r(i,4)*180/3.1415926535;
%     axis([-maximum maximum -maximum maximum]);
% end

%Length shadow

% for i = [1:l(1)]
%     p = r(i,1)*r(i,1)+r(i,2)*r(i,2)+r(i,3)*r(i,3);
%     if (p>20)
%         p = 20;
%     end
%     
%     %plot(r(i,3),r(i,1),'r.');
%     
%     h = r(i,5)/3600.0;
%     plot(h,p,'r.');
%     %axis([-maximum maximum -maximum maximum]);
%     pause(0.001);
% end

%Hist daylight

r = read_data("res2.txt",2);

l = size(r);
day = [1:365];

subplot(2,2,1);
hold on
plot(day,r(:,1)./3600.0,'-r.');
plot(day,r(:,2)./3600.0,'-r.');
plot(day, ones(1,365)*8, '-b.');
plot(day, ones(1,365)*20, '-b.');
axis([1 365 0 24])

subplot(2,2,2);
hold on
plot(day,r(:,1)./3600.0+1,'-r.');
plot(day,r(:,2)./3600.0+1,'-r.');
plot(day, ones(1,365)*8, '-b.');
plot(day, ones(1,365)*20, '-b.');
axis([1 365 0 24])

winter_day = [1:85, 301:365];
summer_day = [86:300];

subplot(2,2,3);
hold on
plot(winter_day,r(winter_day,1)./3600.0,'r.');
plot(winter_day,r(winter_day,2)./3600.0,'r.');

plot(summer_day,r(summer_day,1)./3600.0+1,'r.');
plot(summer_day,r(summer_day,2)./3600.0+1,'r.');
plot(day, ones(1,365)*8, '-b.');
plot(day, ones(1,365)*20, '-b.');
axis([1 365 0 24])

min_h = 8;
max_h = 20;

subplot(2,2,4)
hold on
plot(day,min(r(:,2)./3600.0,max_h) - max(r(:,1)./3600.0,min_h),'-r.');
plot(day,min(r(:,2)./3600.0+1,max_h) - max(r(:,1)./3600.0+1,min_h),'-y.');
plot(winter_day,min(r(winter_day,2)./3600.0,max_h) - max(r(winter_day,1)./3600.0,min_h),'-b.');
plot(summer_day,min(r(summer_day,2)./3600.0+1,max_h) - max(r(summer_day,1)./3600.0+1,min_h),'-b.');
axis([1 365 0 24])