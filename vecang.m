%% Calculates minimum angle between two vectors

function theta=vecang(x,y)

theta = acosd( dot(x,y)/norm(x)/norm(y) );

theta=real(theta);