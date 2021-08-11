function ifdo(condition, eval_if_true, eval_if_false)

if(condition)
    eval_if_true();
else
    if(exist('eval_if_false','var') && ~isempty(eval_if_false))
        eval_if_false();
    end
end
