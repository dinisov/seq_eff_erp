function [PHOT] = trim_phot_outliers(PHOT, photSD)
    %     trim horrible outliers from photodiode data

    PHOT = normalize(PHOT.').';

    PHOT(1,PHOT(1,:) > photSD) = photSD;
    PHOT(1,PHOT(1,:) < -photSD) = -photSD;
    PHOT(2,PHOT(2,:) > photSD) = photSD;
    PHOT(2,PHOT(2,:) < -photSD) = -photSD;
    if size(PHOT,1) > 2
        PHOT(3,PHOT(3,:) > photSD) = photSD;
        PHOT(3,PHOT(3,:) < -photSD) = -photSD;
    end
end

