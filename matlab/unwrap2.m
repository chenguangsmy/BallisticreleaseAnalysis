function val=unwrap2(angle,threshold,mode)
    if(strcmp(mode,'up'))
        if(angle<threshold)
            val=angle+2*pi;
        else
            val=angle;
        end
    elseif(strcmp(mode,'down'))
        if(angle>threshold)
            val=angle-2*pi;
        else
            val=angle;
        end
    end
end
