function m = quickmedian(M) %#codegen

low = 1;
high = size(M,1);

medianidx = fix((low + high) / 2);


%when case 1 works, low
arr = M;
while(1)
    
    %/* One element only */
    if (high <= low)
        break;
    end
    
    %/* Two elements only */
    if (high == low + 1)
        if (arr(low) > arr(high))
            [arr(low) arr(high)] = swapaux(arr(low), arr(high));
        end
        break;
    end
    
    %/* Find medianidx of low, middle and high items; swap to low position */
    middle = fix((low+high)/2);
    
    if (arr(middle) > arr(high))
        [arr(middle) arr(high)] = swapaux(arr(middle), arr(high));
    end
    if (arr(low) > arr(high))
        [arr(low) arr(high)] = swapaux(arr(low), arr(high));
    end
    if (arr(middle) > arr(low))
        [arr(middle) arr(low)] = swapaux(arr(middle), arr(low));
    end
    
    %/* Swap low item (now in position middle) into position (low+1) */
    [arr(middle) arr(low+1)] = swapaux(arr(middle), arr(low+1));
    
    %/* Work from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    while(1)
        ll = ll + 1;
        while(arr(low) > arr(ll))
            ll = ll + 1;
        end
        hh = hh - 1;
        while(arr(hh) > arr(low))
            hh = hh - 1;
        end
        
        if (hh < ll)
            break;
        end
        [arr(ll) arr(hh)] = swapaux(arr(ll), arr(hh));
    end
    
    %/* Swap middle item (in position low) back into correct position */
    [arr(low) arr(hh)] = swapaux(arr(low), arr(hh));
    
    %/* Reset active partition */
    if (hh <= medianidx)
        low = ll;
    end
    if (hh >= medianidx)
        high = hh - 1;
    end
end

m = arr(medianidx);

function [b a] = swapaux(a,b) %#codegen
