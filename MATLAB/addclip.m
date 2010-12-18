
function imorigin = addclip(Xlogo, Ylogo, ...
    feat_logo, imorigin, imclip, imtest, bw, mask)
% imorigin : RGB
% imclip: RGB
% imtest: grayscale
% imlogo: grayscale
% bw: binary
corner_test = corner_detect(imtest);

[Ytest, Xtest] = vanms( corner_test, 300);
Ytest = max(1, floor(Ytest));
Xtest = max(1, floor(Xtest));

feat_test = HOG( imtest, Ytest, Xtest);

match = feat_match( feat_test, feat_logo);
%for each test pt find logo match

index = find( match >= 0);
index_logo = match(index);

m_test_pts = [Xtest(index), Ytest(index)];
m_logo_pts = [Xlogo(index_logo), Ylogo(index_logo)];

[model inlier] = ransac_est_model( m_test_pts(:,1), m_test_pts(:,2), ...
    m_logo_pts(:,1), m_logo_pts(:,2), 1, 1000);
% transform pts on test to logo template and retrive pixel value back

rows = ind2sub( [size(m_test_pts,1), 1], inlier );
in_test_pts = m_test_pts(rows, :);

imorigin = sample_from_clipart(imorigin, in_test_pts, ...
                        bw, imclip, model, mask);

