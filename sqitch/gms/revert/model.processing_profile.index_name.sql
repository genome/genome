-- Revert model.processing_profile.index_name

BEGIN;

DROP INDEX model.processing_profile_name_index;

COMMIT;
