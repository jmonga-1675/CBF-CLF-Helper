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
        %config
               
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
        % params: dictionary of necessary parameters to define the dynsys
        %   you can pass the model parameters you want to use to define
        %   the dynamics, together with the below fields that is necessary
        %   to support the functionality of the library.
        % Necessary fields:
        %   clf_rate / clf.rate: rate for the clf constraint.
        %   cbf_rate / cbf.rate: rate for the cbf constraint.
        % for the 'built-in' option:
        %   xdim: state dimension
        %   udim: control input dimension
        % Optional fields:
        %   u_min: control min bound (scalar, or the vector of size u_min)
        %   u_max: control max bound        
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
        %   x can be multiple element (size: (xdim, n_element))
        %   u can be multiple element (size: (udim, n_element))       
        % Output: dx: \dot(x) (size: (xdim, n_element))
            dx = obj.f(x) + obj.g(x) * u;
        end
        
        function f_ = f(obj, x)
        % f_ = obj.f(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of f(x).
        % Autonomous vector fields (or drift).
        % x can be multiple element (size: (xdim, n_element))
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, f(x) should be overriden by user.");                
            end
            f_ = obj.f_sym(x);            
        end
        
        function g_ = g(obj, x)
        % g_ = obj.g(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of g(x).
        % Control vector fields (or actuation effect).
        % x can be multiple element (size: (xdim, n_element))
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.g(x) should be overriden by user.");                
            end
            g_ = obj.g_sym(x);
        end
        
        function Vs = clf(obj, x)
        % Vs = obj.clf(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the Control Lyapunov Function V(x).
        % x can be multiple elements (size: (xdim, n_element))
        % Vs:  (size: (obj.n_clf, n_element))
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.clf(x) should be overriden by user.");                
            end
            n_states = size(x, 2);
            Vs = zeros(obj.n_clf, n_states);
            for i = 1:obj.n_clf
                clf_i = obj.clf_sym{i}(x);
                Vs(i, :) = clf_i;
            end
        end
        
        function LfVs = lf_clf(obj, x)
        % LfVs = obj.lf_clf(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CLF L_f{V(x)}.
        % x can be multiple elements (size: (xdim, n_element))
        % LfVs:  (size: (obj.n_clf, n_element))
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lf_clf(x) should be overriden by user.");
            end
            n_states = size(x, 2);
            LfVs = zeros(obj.n_clf, n_states);
            for i = 1:obj.n_clf               
                lf_clf_i = obj.lf_clf_sym{i}(x);
                LfVs(i, :) = lf_clf_i;
            end
        end
        
        function LgVs = lg_clf(obj, x)
        % LgVs = obj.lg_clf(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CLF L_g{V(x)}.
        % x can be multiple elements (size: (xdim, n_element))
        % LgVs:  (size: (obj.n_clf, obj.udim, n_element))        
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lg_clf(x) should be overriden by user.");
            end
            n_states = size(x, 2);
            LgVs = zeros(obj.n_clf, obj.udim, n_states);
            for i = 1:obj.n_clf
                lg_clf_i = obj.lg_clf_sym{i}(x);
                LgVs(i, :, :) = lg_clf_i;
            end
        end
        
        function Vdots = dclf(obj, x, u)
        % Vdots = obj.dclf(x, u)
        % Model based estimate of the Lie derivatives of CLF
        % x can be multiple elements (size: (xdim, n_element))
        % u can be multiple elements (size: (udim, n_element))
        % Output size: (obj.n_clf, n_element)
            u = reshape(u, [1, size(u)]);
            u = permute(u, [2, 1, 3]);
            LfVs = obj.lf_clf(x);
            LgVs = obj.lg_clf(x);
            Vdots = LfVs + reshape(pagemtimes(LgVs, u), obj.n_clf, []);            
        end
            
        function Bs = cbf(obj, x)
        % Bs = obj.cbf(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the Control Barrier Function B(x).
        % x can be multiple elements (size: (xdim, n_element))
        % Bs:  (size: (obj.n_cbf, n_element))
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.cbf(x) should be overriden by user.");                
            end
            n_states = size(x, 2);
            Bs = zeros(obj.n_cbf, n_states);
            for i = 1:obj.n_cbf
                cbf_i = obj.cbf_sym{i}(x);
                Bs(i, :) = cbf_i;
            end
        end
        
        function LfBs = lf_cbf(obj, x)
        % LfBs = obj.lf_cbf(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CBF L_f{B(x)}.
        % x can be multiple elements (size: (xdim, n_element))
        % LfBs: (size: (obj.n_cbf, n_element))
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lf_cbf(x) should be overriden by user.");
            end
            n_states = size(x, 2);
            LfBs = zeros(obj.n_cbf, n_states);
            for i = 1:obj.n_cbf               
                lf_cbf_i = obj.lf_cbf_sym{i}(x);
                LfBs(i, :) = lf_cbf_i;
            end
        end
        
        function LgBs = lg_cbf(obj, x)
        % LgBs = obj.lg_cbf(x)
        % For 'built-in' setup, override this function with the
        % user-defined implementation of the lie derivative of the CBF L_g{B(x)}.
        % x can be multiple elements (size: (xdim, n_element))
        % LgBs: (size: (obj.n_cbf, obj.udim, n_element))
            if strcmp(obj.setup_option, 'built-in')
                error("For 'built-in' setup_option, obj.lg_cbf(x) should be overriden by user.");
            end
            n_states = size(x, 2);
            LgBs = zeros(obj.n_cbf, obj.udim, n_states);
            for i = 1:obj.n_cbf
                lg_cbf_i = obj.lg_cbf_sym{i}(x);
                LgBs(i, :, :) = lg_cbf_i;
            end
        end
        
        function Bdots = dcbf(obj, x, u)
        % Bdots = obj.dcbf(x, u)
        % Model based estimate of the Lie derivatives of CBF
        % x can be multiple elements (size: (xdim, n_element))
        % u can be multiple elements (size: (udim, n_element))
        % Output size: (obj.n_cbf, n_element)
            u = reshape(u, [1, size(u)]);
            u = permute(u, [2, 1, 3]);
            LfBs = obj.lf_cbf(x);
            LgBs = obj.lg_cbf(x);
            Bdots = LfBs + reshape(pagemtimes(LgBs, u), obj.n_cbf, []);            
        end
    end
end

