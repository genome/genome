-- Deploy config.analysis_project.index_analysis_menu_item_id
-- requires: analysis_project

BEGIN;

CREATE INDEX c_ap_analysis_menu_item_id_index ON config.analysis_project USING btree (analysis_menu_item_id);

COMMIT;
