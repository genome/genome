-- Revert config.analysismenu_item_indices

BEGIN;
DROP INDEX config.analysismenu_item_name_idx;
DROP INDEX config.analysismenu_item_path_idx;
COMMIT;
