# Here we define some kernels for Smoluchowski


# -------- Smoluchowski: --------------
function ω_Smoluchowski(k)
    return k
end

function K_Smoluchowski_one(k, p)
    return 1
end

function K_Smoluchowski_additive(k, p)
    return k + p
end

function K_Smoluchowski_multiplicative(k, p)
    return k * p
end
