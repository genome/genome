-- Revert result.user.index_user_id_user_class_name_label

BEGIN;

DROP INDEX result.user_id_name_label_index;

COMMIT;
