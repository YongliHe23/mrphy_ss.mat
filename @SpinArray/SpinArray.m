% object for MRI simulations.
%
% Brief intro on methods:
%
%   Tianrui Luo, 2017
%}

classdef SpinArray < matlab.mixin.SetGet & matlab.mixin.Copyable
  methods (Static) % for quickly retriving attributes names. zjdsfyzs
    function n = immutable(), n={'dim', 'mask', 'nM'}; end
    function n = compact(),   n={'M_', 'T1_', 'T2_', 'gam_'}; end
    function n = dependent(), n={'M', 'T1', 'T2', 'gam'}; end
  end

  properties (SetAccess = immutable)
    dim;  % [nx, [ny, [nz, ...]]], i.e., Nd
    mask; % (*Nd) a.u., simulation mask,
    nM;
  end
  properties (SetAccess = public, GetAccess = public)
    M_   % (nM, xyz) spin vector, rotating frame
    T1_  % (nM,) Sec, dflt for grey matter
    T2_  % (nM,) Sec, dflt for grey matter
    gam_ % (nM,) Hz/G
  end
  properties (Dependent = true)
    M
    T1
    T2
    gam
  end

  methods (Access = public)
    function obj = SpinArray(dim, varargin)
      %INPUTS:
      % - dim (nx, (ny, (nz, ...)))
      %OPTIONALS
      % - mask (*Nd) logical
      % - M (*Nd, xyz)
      % - T1 (1,), global; (*Nd, ), spin-wise
      % - T2 (1,), global; (*Nd, ), spin-wise
      % - gam (1,), global; (*Nd, ), spin-wise
      import attr.*

      obj.dim = dim;

      %% dflts
      kv_c = [{'mask'}, obj.dependent; cell(1, 1+numel(obj.dependent))];
      arg = struct(kv_c{:});
      arg.mask = true(obj.dim);
      arg.T1  = ones(obj.dim)*1.47;
      arg.T2  = ones(obj.dim)*0.07;
      arg.gam = ones(obj.dim)*mrphy.utils.envMR('get','gam');
      arg.M   = cat(numel(obj.dim)+1, zeros([obj.dim, 2]), ones(obj.dim));

      arg = attrParser(arg, varargin);

      %%
      [obj.mask, obj.nM] = deal(arg.mask, nnz(arg.mask));
      obj.T1_  = obj.extract(arg.T1);
      obj.T2_  = obj.extract(arg.T2);
      obj.gam_ = obj.extract(arg.gam);
      obj.M_   = obj.extract(arg.M);
    end

    function [Mo_, Mhst_] = applypulse(obj, pulse, varargin)
      %INPUTS:
      % - pulse (1,) @Pulse
      %OPTIONALS:
      % - loc_   ^ loc   (nM, xyz)       ^ (*Nd, xyz), XOR
      % - b0Map_ | b0Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
      % - b1Map_ | b1Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
      % - doCim [T/f]
      % - doEmbed [t/F]
      % - doUpdate [t/F]
      %OUTPUTS:
      % - Mo_   (nM, xyz)     | (*Nd, xyz)
      % - Mhst_ (nM, xyz, nT) | (*Nd, xyz, nT), evolving history
      import attr.*

      %% parsing
      [arg.loc,   arg.loc_] = deal([], []);
      [arg.b0Map, arg.b0Map_] = deal([], []);
      [arg.b1Map, arg.b1Map_] = deal([], []);
      [arg.doCim, arg.doEmbed, arg.doUpdate] = deal(true, false, false);

      arg = attrParser(arg, varargin);

      kw = {  'loc',arg.loc,     'loc_',arg.loc_ ...
            , 'b0Map',arg.b0Map, 'b0Map_',arg.b0Map_ ...
            , 'b1Map',arg.b1Map, 'b1Map_',arg.b1Map_ ...
            , 'doEmbed',false};
      beff_ = obj.pulse2beff(pulse, kw{:});

      [Mo_, Mhst_] = mrphy.sims.blochsim(obj.M_, beff_, obj.T1_, obj.T2_ ...
                                         , pulse.dt, obj.gam_ ...
                                         , arg.doCim);

      if arg.doUpdate, obj.M_ = Mo_; end
      if arg.doEmbed
        [Mo_, Mhst_] = deal(obj.embed(Mo_), obj.embed(Mhst_));
      end
    end
    
    function [M_ss, Mhst_] = applypulse_ss(obj, pulse,varargin)
      %Comupte steady state Magnetization
      %INPUTS:
      % - pulse (1,) @Pulse
      %OPTIONALS:
      % - loc_   ^ loc   (nM, xyz)       ^ (*Nd, xyz), XOR
      % - b0Map_ | b0Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
      % - b1Map_ | b1Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
      % - doCim [T/f]
      % - doEmbed [t/F]
      % - doUpdate [t/F]
      %OUTPUTS:
      % - Mo_   (nM, xyz)     | (*Nd, xyz)
      % - Mhst_ (nM, xyz, nT) | (*Nd, xyz, nT), evolving history
      import attr.*

      %% parsing
      [arg.loc,   arg.loc_] = deal([], []);
      [arg.b0Map, arg.b0Map_] = deal([], []);
      [arg.b1Map, arg.b1Map_] = deal([], []);
      [arg.doCim, arg.doEmbed, arg.doUpdate] = deal(true, false, false);
      [arg.alpha,arg.TR] = deal(15,55e-3); %\degree, ms

      arg = attrParser(arg, varargin);

      kw = {  'loc',arg.loc,     'loc_',arg.loc_ ...
            , 'b0Map',arg.b0Map, 'b0Map_',arg.b0Map_ ...
            , 'b1Map',arg.b1Map, 'b1Map_',arg.b1Map_ ...
            , 'doEmbed',false};
      beff_ = obj.pulse2beff(pulse, kw{:});
      
      %[Mo_, Mhst_] = mrphy.sims.blochsim(obj.M, beff_, obj.T1_, obj.T2_ ...
      [Mo_, Mhst_] = mrphy.sims.blochsim(obj.M_, beff_, obj.T1_, obj.T2_ ...
                                         , pulse.dt, obj.gam_ ...
                                         , arg.doCim);

      alpha = arg.alpha;
      E1=exp(-arg.TR./obj.T1_);

      Mo_(isnan(Mo_))=0; %nan to 0
      beta_=atan((Mo_(:,1).^2+Mo_(:,2).^2).^0.5./Mo_(:,3));

      M_ss=zeros(size(Mo_));

      M_ss(:,3)=obj.M_(:,3).*(1-E1)./(1-cos(beta_).*cosd(alpha).*E1).*cos(beta_);
      M_ss(:,1)=obj.M_(:,3).*(1-E1)./(1-cos(beta_).*cosd(alpha).*E1).*sin(beta_);

      if arg.doUpdate, obj.M_ = M_ss; end
      if arg.doEmbed
        [M_ss, Mhst_] = deal(obj.embed(M_ss), obj.embed(Mhst_));
      end
    end

    function [M_ss, Mhst_] = applypulse_ss_sms(obj, pulse, varargin)
      %Compute steady state Magnetization for inner-volume saturation SS-SMS EPI
      %Sequence: [beta(iv_saturate)--alpha(slice-selective_tip_down)-readout]
      %          + [beta--readout] x (vTR/TR - 1)
      %INPUTS:
      % - pulse (1,) @Pulse
      %OPTIONALS:
      % - loc_   ^ loc   (nM, xyz)       ^ (*Nd, xyz), XOR
      % - b0Map_ | b0Map (nM, 1)         ^ (*Nd, 1)
      % - b1Map_ | b1Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
      % - doCim  [T/f]
      % - doEmbed [t/F]
      % - doUpdate [t/F]
      % - TR (1,) s, repetition time for beta pulse. Dflt 55e-3.
      % - vTR (1,) s, volume TR (repetition time for alpha). Dflt 55e-2.
      % - alpha (1,) deg, flip angle of slice-selective tip-down. Dflt 52.
      % - alphaDur (1,) s, duration of alpha pulse. Dflt 8e-3.
      %OUTPUTS:
      % - M_ss  (nM, xyz) | (*Nd, xyz)
      % - Mhst_ (nM, xyz, nT) | (*Nd, xyz, nT)
      import attr.*

      %% parsing
      [arg.loc,   arg.loc_] = deal([], []);
      [arg.b0Map, arg.b0Map_] = deal([], []);
      [arg.b1Map, arg.b1Map_] = deal([], []);
      [arg.doCim, arg.doEmbed, arg.doUpdate] = deal(true, false, false);
      [arg.alpha, arg.TR, arg.vTR, arg.alphaDur] = deal(52, 55e-3, 55e-2, 8e-3);

      arg = attrParser(arg, varargin);

      kw = {  'loc',arg.loc,     'loc_',arg.loc_ ...
            , 'b0Map',arg.b0Map, 'b0Map_',arg.b0Map_ ...
            , 'b1Map',arg.b1Map, 'b1Map_',arg.b1Map_ ...
            , 'doEmbed',false};
      [Mo_, Mhst_] = obj.applypulse(pulse, kw{:});

      nT = size(pulse.rf, 2);
      Tr = arg.TR - pulse.dt * nT - arg.alphaDur; % time after beta+alpha to end of TR
      Tt = arg.TR - pulse.dt * nT;               % time after beta to end of TR

      Mo_(isnan(Mo_)) = 0;
      Beta_ = atan(sqrt(Mo_(:,1).^2 + Mo_(:,2).^2) ./ Mo_(:,3)); % (nM,) rad

      E1_  = exp(-Tr ./ obj.T1_); % T1 relaxation over Tr
      E1t_ = exp(-Tt ./ obj.T1_); % T1 relaxation over Tt

      A_  = cos(Beta_) .* E1t_;
      B_  = obj.M_(:,3) .* (1 - E1t_);
      As_ = cos(Beta_) .* E1_ .* cosd(arg.alpha); % A*
      Bs_ = obj.M_(:,3) .* (1 - E1_);             % B*

      Nrep  = arg.vTR / arg.TR;
      Nrep_ = Nrep - 1;

      M_ss = zeros(size(Mo_));
      M_ss(:,3) = (A_.^Nrep_ .* Bs_ + (1 - A_.^Nrep_) .* B_ ./ (1 - A_)) ...
                  ./ (1 - A_.^Nrep_ .* As_);

      if arg.doUpdate, obj.M_ = M_ss; end
      if arg.doEmbed
        [M_ss, Mhst_] = deal(obj.embed(M_ss), obj.embed(Mhst_));
      end
    end

    function st = asstruct(obj)
      warning('off', 'MATLAB:structOnObject')
      st = struct(obj);
      warning('on', 'MATLAB:structOnObject')
    end

    function target = target_Mss(obj, beta_iv, beta_ov, iv, ov, varargin)
      %Compute target steady-state magnetization profile before alpha excitation
      %INPUTS:
      % - beta_iv (1,) deg, target iv beta flip angle
      % - beta_ov (1,) deg, target ov beta flip angle
      % - iv (*Nd) logical, inner volume mask
      % - ov (*Nd) logical, outer volume mask
      %OPTIONALS:
      % - weight_iv (1,) [1.0]
      % - weight_ov (1,) [1.0]
      % - doEmbed [T/f]
      % - alpha   (1,) deg [15]
      % - TR      (1,) Sec [55e-3]
      %OUTPUTS:
      % - target struct with fields:
      %     d      (nM, xyz) | (*Nd, xyz), target SS magnetization
      %     weight (nM,)     | (*Nd,),    voxel weights
      import attr.*

      [arg.weight_iv, arg.weight_ov] = deal(1.0, 1.0);
      [arg.doEmbed, arg.alpha, arg.TR] = deal(true, 15, 55e-3);
      arg = attrParser(arg, varargin);

      iv_ = obj.extract(iv);  % (nM, 1)
      ov_ = obj.extract(ov);  % (nM, 1)

      M0 = obj.M_(:, 3);      % (nM, 1), z-component
      E1 = exp(-arg.TR ./ obj.T1_);  % (nM, 1)

      beta_rad  = deg2rad(beta_iv .* iv_ + beta_ov .* ov_);  % (nM, 1)
      alpha_rad = deg2rad(arg.alpha);

      d_ = zeros(obj.nM, 3);
      d_(:, 3) = (M0 .* (1 - E1)) ./ (1 - cos(beta_rad) .* cos(alpha_rad) .* E1) .* cos(beta_rad);
      d_(:, 1) = (M0 .* (1 - E1)) ./ (1 - cos(beta_rad) .* cos(alpha_rad) .* E1) .* sin(beta_rad);
      d_(isnan(d_)) = 0;

      weight_ = arg.weight_iv .* ov_ + arg.weight_ov .* iv_;  % (nM, 1)

      if arg.doEmbed
        target = struct('d', obj.embed(d_), 'weight', obj.embed(weight_));
      else
        target = struct('d', d_, 'weight', weight_);
      end
    end

    function Mr_ = freeprec(obj, dur, varargin)
      %INPUTS:
      % - dur (1,) Sec, duration of free precision
      %OPTIONALS:
      % - b0Map_ | b0Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils), "Hz"
      % - doEmbed [t/F]
      % - doUpdate [t/F]
      %OUTPUTS:
      % - Mr_ (nM, xyz) | (*Nd, xyz)
      import attr.*

      %% parsing
      [arg.b0Map, arg.b0Map_] = deal([], []);
      [arg.doEmbed, arg.doUpdate] = deal(false, false);

      arg = attrParser(arg, varargin);

      assert( isempty(arg.b0Map) || isempty(arg.b0Map_) ) % Not both
      if ~isempty(arg.b0Map), arg.b0Map_ = obj.extract(arg.b0Map); end

      Mr_ = freePrec(obj.M_, dur, obj.T1_, obj.T2_, arg.b0Map_);

      if arg.doEmbed, Mr_ = obj.embed(Mr_); end
      if arg.doUpdate, obj.M_ = Mr_; end
    end

    function beff_ = pulse2beff(obj, pulse, varargin)
      %INPUTS:
      % - pulse (1,) mrphy.@Pulse
      %OPTIONALS:
      % - loc_   ^ loc   (nM, xyz)       ^ (*Nd, xyz), XOR
      % - b0Map_ | b0Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
      % - b1Map_ | b1Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
      % - doEmbed [t/F]
      %OUTPUTS:
      % - beff_ (nM, xyz, nT) | (*Nd, xyz, nT)
      import attr.*

      %% parsing
      [arg.loc,   arg.loc_] = deal([], []);
      [arg.b0Map, arg.b0Map_] = deal([], []);
      [arg.b1Map, arg.b1Map_] = deal([], []);
      arg.doEmbed = false;

      arg = attrParser(arg, varargin);

      assert( isempty(arg.loc) ~= isempty(arg.loc_) ) % XOR
      if isempty(arg.loc_), arg.loc_ = obj.extract(arg.loc); end

      assert( isempty(arg.b0Map) || isempty(arg.b0Map_) ) % Not both
      if ~isempty(arg.b0Map), arg.b0Map_ = obj.extract(arg.b0Map); end

      assert( isempty(arg.b1Map) || isempty(arg.b1Map_) ) % Not both
      if ~isempty(arg.b1Map), arg.b1Map_ = obj.extract(arg.b1Map); end

      %%
      kw = {'gam', obj.gam_, 'b0Map', arg.b0Map_, 'b1Map', arg.b1Map_};
      beff_ = pulse.beff(arg.loc_, kw{:});
      if arg.doEmbed, beff_ = obj.embed(beff_); end
    end

  end
  
  methods % Utilities
    function v_ = extract(obj, v)
      [mask_t, ndim] = deal(obj.mask, numel(obj.dim));
      s_v = size(v);
      if numel(s_v) < ndim, s_v = [s_v, ones(1, ndim-numel(s_v))]; end
      shape_v = [s_v, 1]; % [*Nd, ..., 1]
      v = reshape(v, prod(obj.dim), []);
      v_ = reshape(v(mask_t,:),[obj.nM,shape_v(ndim+1:end)]); % (nM, ...)
    end

    function v = embed(obj, v_)
      shape_v_ = [size(v_), 1];
      v = nan([prod(obj.dim), shape_v_(2:end)]);
      v(obj.mask,:) = v_(:,:);
      v = reshape(v, [obj.dim, shape_v_(2:end)]);
    end
  end

  methods % set and get, sealed if the property cannot be redefined
    %% Dependent variables
    % WARNING
    % DO NOT proceed indexed/masked assignment to non-compact properties.
    function set.T1(obj, v),  obj.T1_  = obj.extract(v); end
    function set.T2(obj, v),  obj.T2_  = obj.extract(v); end
    function set.gam(obj, v), obj.gam_ = obj.extract(v); end
    function set.M(obj, v),   obj.M_   = obj.extract(v); end

    function v = get.T1(obj),  v = obj.embed(obj.T1_);  end
    function v = get.T2(obj),  v = obj.embed(obj.T2_);  end
    function v = get.gam(obj), v = obj.embed(obj.gam_); end
    function v = get.M(obj),   v = obj.embed(obj.M_);   end
  end
end
