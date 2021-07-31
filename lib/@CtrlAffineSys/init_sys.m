function init_sys(obj, params)
    obj.params = params;

    if strcmp(obj.setup_option, 'symbolic')
        disp(['Setting up the dynamics, CLFs, CBFs from defined symbolic expressions.', ...
            '(This might take time.)']);
        % get symbolic expressions for x, f, g, clf, cbf.
        [x, f_, g_] = obj.defineSystem(params);
        if isempty(x) || isempty(f_) || isempty(g_)
            error("x, f, g is empty. Create a class function defineSystem %s", ...
            "and define your dynamics with symbolic expression.");
        end
        if ~isa(f_, 'sym')
            f_ = sym(f_);
        end
        if ~isa(g_, 'sym')
            g_ = sym(g_);
        end
        clf_ = obj.defineClf(params, x);
        cbf_ = obj.defineCbf(params, x);
        % Setting state and input dimension.
        obj.xdim = size(x, 1);
        obj.udim = size(g_, 2);
        obj.f_sym = matlabFunction(f_, 'vars', {x});
        obj.g_sym = matlabFunction(g_, 'vars', {x});
        
        if ~isempty(clf_)
            if ~iscell(clf_)
                if length(clf_) ~=1
                    error("Defined CLF should be a scalar. For multiple CLFs, %s", ...
                        "use a cell array.");
                end
                clf_ = {clf_};
            end
            obj.n_clf = length(clf_);
            clf_sym = cell(obj.n_clf, 1);
            lf_clf_sym = cell(obj.n_clf, 1);
            lg_clf_sym = cell(obj.n_clf, 1);
            for i_clf = 1:obj.n_clf
                dclf_i = simplify(jacobian(clf_{i_clf}, x));
                lf_clf_i = dclf_i * f_;
                lg_clf_i = dclf_i * g_;
                clf_sym{i_clf} = matlabFunction(clf_{i_clf}, 'vars', {x});
                lf_clf_sym{i_clf} = matlabFunction(lf_clf_i, 'vars', {x});
                if all(isAlways(simplify(lg_clf_i) == 0, 'Unknown', 'false'))
                    error('Relative degree of the defined %d-th CLF > 1. %s', ...
                    [i_clf, 'Currently, High-order relative degree is not supported.']);
                end
                lg_clf_sym{i_clf} = matlabFunction(lg_clf_i, 'vars', {x});
            end
            obj.clf_sym = clf_sym;             
            obj.lf_clf_sym = lf_clf_sym;
            obj.lg_clf_sym = lg_clf_sym;
        else
            obj.n_clf = 0;
        end
        if ~isempty(cbf_)
            if ~iscell(cbf_)
                if length(cbf_) ~=1
                    error("Defined CBF should be a scalar. For multiple CBFs, %s", ...
                        "use a cell array.");
                end
                cbf_ = {cbf_};
            end
            obj.n_cbf = length(cbf_);
            cbf_sym = cell(obj.n_cbf, 1);
            lf_cbf_sym = cell(obj.n_cbf, 1);
            lg_cbf_sym = cell(obj.n_cbf, 1);
            for i_cbf = 1:obj.n_cbf
                dcbf_i = simplify(jacobian(cbf_{i_cbf}, x));
                lf_cbf_i = dcbf_i * f_;
                lg_cbf_i = dcbf_i * g_;        
                cbf_sym{i_cbf} = matlabFunction(cbf_{i_cbf}, 'vars', {x});
                lf_cbf_sym{i_cbf} = matlabFunction(lf_cbf_i, 'vars', {x});
                if all(isAlways(simplify(lg_cbf_i) == 0, 'Unknown', 'false'))
                    error('Relative degree of the defined %d-th CBF > 1. %s', ...
                    [i_cbf, 'Currently, High-order relative degree is not supported.']);
                end
                lg_cbf_sym{i_cbf} = matlabFunction(lg_cbf_i, 'vars', {x});                
            end
            obj.cbf_sym = cbf_sym;             
            obj.lf_cbf_sym = lf_cbf_sym;
            obj.lg_cbf_sym = lg_cbf_sym;
        else
            obj.n_cbf = 0;
        end
        
    elseif strcmp(obj.setup_option, 'built-in')
        if ~isfield(params, 'xdim')
            error("xdim should be specified for built-in setup.");
        end
        obj.xdim = params.xdim;
        if ~isfield(params, 'udim')
            error("udim should be specified for built-in setup.");
        end
        obj.udim = params.udim;
        x_test = zeros(obj.xdim, 1);
        n_clf = 1;
        try
            obj.clf(x_test);
        catch e
            n_clf = 0;
        end
        obj.n_clf = n_clf;
        n_cbf = 1;
        try
            obj.cbf(x_test);
        catch e
            n_cbf = 0;
        end
        obj.n_cbf = n_cbf;               
    else
        error("Undefined setup_option.");
    end
    
    %% Parse parameters for both setup_option.
    if isfield(params, 'clf') && isfield(params.clf, 'rate')
        %% TODO: change this to params.clf_rate
        obj.clf_rate = params.clf.rate;
    elseif isfield(params, 'clf_rate')
        obj.clf_rate = params.clf_rate;
    elseif obj.n_clf >= 1
        error("params.clf.rate or params.clf_rate should be provided %s", ...
            "in order to use CLF.");
    end
    if isfield(params, 'cbf') && isfield(params.cbf, 'rate')
        %% TODO: change this to params.cbf_rate            
        obj.cbf_rate = params.cbf.rate;
    elseif isfield(params, 'cbf_rate')
        obj.cbf_rate = params.cbf_rate;
    elseif obj.n_cbf >= 1
        error("params.cbf.rate or params.cbf_rate should be provided %s", ...
            "in order to use CBF.");
    end
    if isfield(params, 'u_min')
        if length(params.u_min) == 1
            obj.u_min = params.u_min * ones(obj.udim, 1);
        elseif length(params.u_min) ~= obj.umin
            error("Invalid size of params.u_min.");
        else
            if isrow(params.u_min)
                obj.u_min = params.u_min';
            else
                obj.u_min = params.u_min;
            end
        end
    end
    if isfield(params, 'u_max')
        if length(params.u_max) == 1
            obj.u_max = params.u_max * ones(obj.udim, 1);
        elseif length(params.u_max) ~= obj.umax
            error("Invalid size of params.u_max.");
        else
            if isrow(params.u_max)
                obj.u_max = params.u_max';
            else
                obj.u_max = params.u_max;
            end
        end
    end
    
    %% Do sanity check if all necessary functions are set up properly.
    x_test = zeros(obj.xdim, 1);
    % Vector fields.
    % Error should occure if they are not set properly.
    f_test = obj.f(x_test);
    g_test = obj.g(x_test);
    if size(f_test, 1) ~= obj.xdim || size(f_test, 2) ~= 1
        error("f has wrong size.");
    elseif size(g_test, 1) ~= obj.xdim || size(g_test, 2) ~= obj.udim
        error("g has wrong size.");
    end

    if obj.n_clf >= 1
        clf_test = obj.clf(x_test);
        if numel(clf_test) ~= obj.n_clf
            error("clf has wrong size.");
        end
        lf_clf_test = obj.lf_clf(x_test);
        if numel(lf_clf_test) ~= obj.n_clf
            error("lf_clf has wrong size.");
        end
        lg_clf_test = obj.lg_clf(x_test);
        if size(lg_clf_test, 1) ~= obj.n_clf || size(lg_clf_test, 2) ~= obj.udim
            error("lg_clf has wrong size.");
        end
    end

    if obj.n_cbf >= 1
        cbf_test = obj.cbf(x_test);
        if numel(cbf_test) ~= obj.n_cbf
            error("cbf has wrong size.");
        end
        lf_cbf_test = obj.lf_cbf(x_test);
        if numel(lf_cbf_test) ~= obj.n_cbf
            error("lf_cbf has wrong size.");
        end
        lg_cbf_test = obj.lg_cbf(x_test);
        if size(lg_cbf_test, 1) ~= obj.n_cbf || size(lg_cbf_test, 2) ~= obj.udim
            error("lg_cbf has wrong size.");
        end
    end    
    
    fprintf(obj.get_dynsys_summary())
end