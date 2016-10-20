function [Xmin,fmin,output] = simple_mads(X0,bb_handle,lb,ub,options)

  %-------------------------------
  % Options
  %-------------------------------

  % Default options
  default_options.budget = +inf;
  default_options.tol    = 1e-9;
  default_options.psize_init = 1.0;
  default_options.display         = false;
  default_options.opportunistic   = true;
  default_options.check_cache     = true;

  % Construct the option structure from user options and default options
  if ~exist('options','var')
      options = struct;
  end
  names = fieldnames(default_options);
  for i = 1:length(names)
      fieldstr = names{i};
      if ~isfield(options,fieldstr)
          options.(fieldstr) = default_options.(fieldstr);
      end
  end

  %-------------------------------
  % Dimension and verifications
  %-------------------------------
  DIM = size(X0,2);
  if length(lb)~=DIM
      error('Wrong lb dimension');
  end
  if length(ub)~=DIM
      error('Wrong ub dimension');
  end

  %-------------------------------
  % Initialization
  %-------------------------------
  iter = 0;
  bbe = 0;
  CACHE = [];
  fmin = +inf;
  hmin = +inf;

  % Variable scaling
  scaling = (ub-lb)/10.0;
  scaling(isinf(scaling)) = 1.0;
  scaling = diag(scaling);

  % Poll size initialization
  psize = options.psize_init;
  psize_success = psize;
  psize_max = 0;

  %-------------------------------
  % Start the optimization
  %-------------------------------
  while true

      if iter==0
          % Evaluate starting points
          if (options.display)
              disp('Evaluation of the starting points');
          end
          POLL = X0;
      else       
          %-------------------------------
          % Build polling directions
          %-------------------------------
          msize = min(psize^2,psize);
          rho = psize/msize;
          % Generate direction
          v = randn(DIM,1);
          % Normalize
          v = v/norm(v);
          % Build Householder matrix
          H = eye(DIM)-2*v*v';
          % Normalization of each column
          H = H*diag(max(abs(H)).^-1);
          % Rounding (and transpose)
          H = msize*ceil(rho*H)';
          % Take the opposite directions
          H = [H;-H];

          % scaling of the directions
          H = H*scaling;

          % Build POLL / central point
          POLL = ones(2*DIM,1)*Xmin + H;

          % Shuffle poll
          POLL = POLL(randperm(2*DIM),:);
      end

      %-------------------------------
      % Evaluate points of the poll
      %-------------------------------
      success = false;
      for i=1:size(POLL,1)

          % Get point
          Xtry = POLL(i,:);

          % Snap to bounds
          Xtry = max(Xtry,lb);
          Xtry = min(Xtry,ub);

          % Search in CACHE
          if options.check_cache && ~isempty(CACHE)
              DC = CACHE-ones(size(CACHE,1),1)*Xtry;
              DC = abs(DC)<options.tol;
              if any(all(DC,2))
                  if options.display
                      disp('Cache hit');
                  end
                  continue;
              end
          end

          % Evaluation
          bbe = bbe+1;
          bbo = bb_handle(Xtry);
          % compute h and f
          [ftry,htry] = eval_fh(bbo);
          
          % Add to the CACHE
          if options.check_cache
              CACHE(end+1,:) = Xtry;
          end

          % Save values
          if ( (hmin>0) && (htry<hmin) ) || ( (htry==0) && (ftry<fmin) )
              success = true;
              fmin = ftry;
              hmin = htry;
              if options.display
                  disp(['Succes : ' num2str(fmin) ' (hmin = ' num2str(hmin) ')']);
              end
              Xmin = Xtry;
              psize_success = psize;
              psize_max = max(psize,psize_max);
          end
          
          if bbe>=options.budget
              break; % Reached the total budget
          end
          if success && options.opportunistic && (iter>1)
              break; % Quit the evaluation of the POLL
          end

      end

      %-------------------------------
      % Updates
      %-------------------------------

      if iter>0
          if success
              psize = psize*2;
          else
              psize = psize/2;
          end
      end

      if options.display
          if iter==0
              disp('End of the evaluation of the starting points');
          end
          disp(['iter=' num2str(iter) ' bbe=' num2str(bbe) ' psize=' num2str(psize,3) '  hmin=' num2str(hmin,3) '  fmin=' num2str(fmin,3)]);
      end

      if (abs(psize)<options.tol) || (bbe>=options.budget)
          break;
      end
      
      iter = iter+1;
  end % end of the mads iterations


  if options.display
      disp(['mads break - iter ' num2str(iter,2) ', psize ' num2str(psize,2) ', bbe ' num2str(bbe)]);
      disp(['Final  : ' num2str(fmin) ' (hmin = ' num2str(hmin) ')']);
  end

  % Build the output
  output.fmin = fmin;
  output.hmin = hmin;
  output.psize = psize;
  output.psize_success = psize_success;
  output.psize_max = psize_max;

end % end function simple_mads



function [f,h] = eval_fh(bbo)
  % Concatenation of bbo from f & c functions
  f = bbo(1);
  if length(bbo)==1
      h = 0;
  else
      c = bbo(2:end);
      % Check non valid values
      if any(isnan(c))
          h = +inf;
      else
          h = sum(max(c,0).^2);
      end
  end
  % Penalize objective
  if isnan(f) || (h>0)
    f = +inf;
  end
end % end function eval_fh
