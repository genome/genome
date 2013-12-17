-- Deploy config.analysismenu_item_indices
-- requires: config_analysismenu_item

BEGIN;
CREATE INDEX analysismenu_item_name_idx ON config.analysismenu_item(name);
CREATE INDEX analysismenu_item_path_idx ON config.analysismenu_item(file_path);
COMMIT;
