function target = target_Mss(cube, beta_iv, beta_ov, iv, ov, varargin)
      %Compute target steady-state magnetization profile before alpha excitation
      %Assuming the following sequence structure:
      %  [beta-alpha-readout]xN_rep
      %where beta is a saturation preparation pulse and alpha is a non-selective excitation pulse.
      %
      %Usage:
      %  ``target = target_Mss(cube, 0, 90, iv_mask, ov_mask, alpha=20, TR=80e-3)``
      %
      %INPUTS:
      % - beta_iv (1,) deg, target inner-volume beta flip angle
      % - beta_ov (1,) deg, target outer-volume beta flip angle
      % - iv (*Nd) logical, inner volume mask
      % - ov (*Nd) logical, outer volume mask
      %OPTIONALS:
      % - weight_iv (1,) [1.0]
      % - weight_ov (1,) [1.0]
      % - doEmbed [T/f]
      % - alpha   (1,) deg [15], uniform excitation pulse FA  
      % - TR      (1,) Sec [55e-3]
      %OUTPUTS:
      % - target struct with fields:
      %     d      (nM, xyz) | (*Nd, xyz), target SS magnetization
      %     weight (nM,)     | (*Nd,),    voxel weights
      import attr.*

      [arg.weight_iv, arg.weight_ov] = deal(1.0, 1.0);
      [arg.doEmbed, arg.alpha, arg.TR] = deal(true, 15, 55e-3);
      arg = attrParser(arg, varargin);

      iv_ = cube.extract(iv);  % (nM, 1)
      ov_ = cube.extract(ov);  % (nM, 1)

      M0 = cube.M_(:, 3);      % (nM, 1), z-component
      E1 = exp(-arg.TR ./ cube.T1_);  % (nM, 1)

      beta_rad  = deg2rad(beta_iv .* iv_ + beta_ov .* ov_);  % (nM, 1)
      alpha_rad = deg2rad(arg.alpha);

      d_ = zeros(cube.nM, 3);
      d_(:, 3) = (M0 .* (1 - E1)) ./ (1 - cos(beta_rad) .* cos(alpha_rad) .* E1) .* cos(beta_rad); %Mz
      d_(:, 1) = (M0 .* (1 - E1)) ./ (1 - cos(beta_rad) .* cos(alpha_rad) .* E1) .* sin(beta_rad); %Mxy (assume along x-axis)
      d_(isnan(d_)) = 0;

      weight_ = arg.weight_iv .* ov_ + arg.weight_ov .* iv_;  % (nM, 1)

      if arg.doEmbed
        target = struct('d', cube.embed(d_), 'weight', cube.embed(weight_));
      else
        target = struct('d', d_, 'weight', weight_);
      end
    end