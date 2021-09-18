%% Author: Jason Choi (jason.choi@berkeley.edu)
classdef CtrlAffineSys < handle    
    %% Control-Affine Dynamic System Class.
    properties
        % Set-up option: 'symbolic', 'built-in'
        setup_option
                
        % necessary parameters
        % rate used in the clf constraint.
        clf_rate
        % rate used in the cbf constraint.
        cbf_rate
        % max bound of control input.
        u_max
        % min bound of control input.
        u_min
                
        % State dimension
        xdim 
        sdim
        % Control input dimension
        udim
        % Number of clf in use
        n_clf
        % Number of cbf in use
        n_cbf
        
        % Model parameters as a structure object, specific to the system.
        % This might contain other parameters used for setting up the dynsys,
        % for instance, clf_rate and cbf_rate. However, referring to this 
        % variable's fields directly instead of the equivalent class variables
        % (when they exist) is not recommended.
        % (Always preferred to refer directly to the class properties.)
        params
        config
        use_phase
               
        %% Functions generated from symbolic expressions.
        % (Used when setup_option is 'symbolic'.)
        f_sym
        g_sym
        cbf_sym
        lf_cbf_sym
        lg_cbf_sym
        clf_sym
        lf_clf_sym
        lg_clf_sym
    end
    
    methods
        function obj = CtrlAffineSys(params, setup_option)
        %% dynsys = CtrlAffineSys(params)
        %% default:
        %%      dynsys = CtrlAffineSys(params, 'symbolic')
        %% if user-defined dynamics, cbf, clf are used:
        %%      dynsys = CtrlAffineSys(params, 'built-in')
            if nargin < 1
                error("Warning: params argument is missing.")
            end
            if nargin < 2
                setup_option = 'symbolic';
            end
            if strcmp(setup_option, 'built_in')
                setup_option = 'built-in';
            elseif strcmp(setup_option, 'builtin')
                setup_option = 'built-in';
            end
            obj.setup_option = setup_option;
            obj.init_sys(params);
        end
        
        function [x, f, g] = defineSystem(obj, params)
        % [x, f, g] = defineSystem(obj, params)
        % For 'symbolic' setup, use this to define your dynamics, 
        % Outputs: x: symbolic state vector
        %          f: drift term, expressed symbolically wrt x.
        %          g: control vector fields, expressed symbolically wrt x.
            x = []; f = []; g = [];         
        end
        
        function clf = defineClf(obj, params, symbolic_state)
        % clf = defineClf(obj, params, symbolic_state)
        % For 'symbolic' setup, use this to define the CLF.
        % symbolic state: same symbolic state created in defineSystem
        % clf: CLF expressed symbolically wrt symbolic_state.
            clf = [];
        end
        
        function cbf = defineCbf(obj, params, symbolic_state)
        % cbf = defineCbf(obj, params, symbolic_state)
        % For 'symbolic' setup, use this to define the CBF.
        % symbolic state: same symbolic state created in defineSystem
        % cbf: CBF expressed symbolically wrt symbolic_state.
            cbf = [];
        end
                
        function dx = dynamics(obj, t, x, u)
        % dx = obj.dynamics(t, x, u)
        % Control-Affine Dynamics: xdot = f(x) + g(x) u
        % Inputs: t: time, x: state, u: control input
        % Output: dx: \dot(x)
            dx = obj.f(x) + obj.g(x) * u;
        end
        
        function f_ = f(obj, x)
        % f_ = f(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of f(x).
        % Autonomous vector fields (or drift).
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, f(x) should be overriden by user.");                
            end
            f_ = obj.f_sym(x);            
        end
        
        function g_ = g(obj, x)
        % g_ = g(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of g(x).
        % Control vector fields (or actuation effect).
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.g(x) should be overriden by user.");                
            end
            g_ = obj.g_sym(x);
        end
        
        function clf_ = clf(obj, x)
        % clf_ = clf(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the Control Lyapunov Function V(x).
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.clf(x) should be overriden by user.");                
            end
            clf_ = zeros(obj.n_clf, 1);
            for i = 1:obj.n_clf
                clf_i = obj.clf_sym{i}(x);
                clf_(i) = clf_i;
            end
        end
        
        function lf_clf_ = lf_clf(obj, x)
        % lf_clf_ = lf_clf(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CLF L_f{V(x)}.
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lf_clf(x) should be overriden by user.");
            end
            lf_clf_ = zeros(obj.n_clf, 1);
            for i = 1:obj.n_clf               
                lf_clf_i = obj.lf_clf_sym{i}(x);
                lf_clf_(i) = lf_clf_i;
            end
        end
        
        function lg_clf_ = lg_clf(obj, x)
        % lg_clf_ = lg_clf(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CLF L_g{V(x)}.
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lg_clf(x) should be overriden by user.");
            end
            lg_clf_ = zeros(obj.n_clf, obj.udim);
            for i = 1:obj.n_clf
                lg_clf_i = obj.lg_clf_sym{i}(x);
                lg_clf_(i, :) = lg_clf_i;
            end
        end
        
        function cbf_ = cbf(obj, x)
        % cbf_ = cbf(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the Control Barrier Function B(x).
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.cbf(x) should be overriden by user.");                
            end
            %cbf_ = zeros(obj.n_cbf, 1);
            cbf_ = [];
            for i = 1:obj.n_cbf
                cbf_i = obj.cbf_sym{i}(x);
                %cbf_(i) = cbf_i;
                cbf_ = [cbf_; cbf_i];
            end
        end
        
        function lf_cbf_ = lf_cbf(obj, x)
        % lf_cbf_ = lf_cbf(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CBF L_f{B(x)}.
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lf_cbf(x) should be overriden by user.");
            end
            lf_cbf_ = zeros(obj.n_cbf, 1);
            for i = 1:obj.n_cbf               
                lf_cbf_i = obj.lf_cbf_sym{i}(x);
                lf_cbf_(i) = lf_cbf_i;
            end
        end
        
        function lg_cbf_ = lg_cbf(obj, x)
        % lg_cbf_ = lg_cbf(obj, x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CBF L_g{B(x)}.
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lg_cbf(x) should be overriden by user.");
            end
            lg_cbf_ = zeros(obj.n_cbf, obj.udim);
            for i = 1:obj.n_cbf
                lg_cbf_i = obj.lg_cbf_sym{i}(x);
                lg_cbf_(i, :) = lg_cbf_i;
            end
        end        
    end
end

