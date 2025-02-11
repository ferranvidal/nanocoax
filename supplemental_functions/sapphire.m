function epsi = sapphire(lambda)


epsi = 1 + 1.4313493*lambda.^2./(lambda.^2-0.0726631^2) + ...
         + 0.65054713*lambda.^2./(lambda.^2-0.1193242^2) + ...
         + 5.3414021*lambda.^2./(lambda.^2-18.028251^2); 