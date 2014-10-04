function [ ] = ecg_display( interval, signal )

plot(interval, signal(1:length(interval)));

end

